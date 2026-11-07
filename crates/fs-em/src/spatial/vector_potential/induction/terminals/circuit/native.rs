//! One native graph preparation and physical account for all electrode fields.
//!
//! The numerical reference remains an interconnection device: its resistor
//! heat cancels against private current-port work and never reaches a physical
//! receiver. Original RLC/transformer/mutual storage stays with CircuitDae.

use super::*;
use crate::spatial::vector_potential::induction::terminals::SpatialInductionTerminalTransfer;

pub(in crate::spatial::vector_potential::induction) struct FieldAccount<'a> {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    pub heat_j: f64,
    pub numerical_dissipation_j: f64,
    pub impressed_work_j: f64,
    pub parameter_work_j: f64,
    pub mechanical_work_j: f64,
    pub terminal_work_j: f64,
    pub terminals: &'a [SpatialInductionTerminalTransfer],
}

/// Only original native field receipts implement this internal view.
pub(in crate::spatial::vector_potential::induction) trait FieldReceipt {
    fn account(&self) -> FieldAccount<'_>;
}
impl FieldReceipt for SpatialInductionStep {
    fn account(&self) -> FieldAccount<'_> {
        FieldAccount {
            stored_before_j: self.stored_before_j,
            stored_after_j: self.stored_after_j,
            stored_change_j: self.stored_after_j - self.stored_before_j,
            heat_j: self.heat_j,
            numerical_dissipation_j: 0.0,
            impressed_work_j: self.impressed_work_j,
            parameter_work_j: self.parameter_work_j,
            mechanical_work_j: self.mechanical_work_j,
            terminal_work_j: self.terminal_work_j,
            terminals: &self.terminals,
        }
    }
}
impl FieldReceipt for super::super::super::nonlinear::SpatialNonlinearInductionStep {
    fn account(&self) -> FieldAccount<'_> {
        FieldAccount {
            stored_before_j: self.stored_before_j,
            stored_after_j: self.stored_after_j,
            stored_change_j: self.stored_change_j,
            heat_j: self.heat_j,
            numerical_dissipation_j: self.numerical_dissipation_j,
            impressed_work_j: self.impressed_work_j,
            parameter_work_j: self.parameter_work_j,
            mechanical_work_j: self.mechanical_work_j,
            terminal_work_j: self.terminal_work_j,
            terminals: &self.terminals,
        }
    }
}
#[cfg(feature = "magnetic-hysteresis")]
impl FieldReceipt for super::super::super::hysteretic::SpatialHystereticInductionStep {
    fn account(&self) -> FieldAccount<'_> {
        FieldAccount {
            stored_before_j: self.stored_before_j,
            stored_after_j: self.stored_after_j,
            stored_change_j: self.stored_change_j,
            heat_j: self.heat_j,
            numerical_dissipation_j: self.numerical_dissipation_j,
            impressed_work_j: self.impressed_work_j,
            parameter_work_j: 0.0,
            mechanical_work_j: 0.0,
            terminal_work_j: self.terminal_work_j,
            terminals: &self.terminals,
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct Allowances {
    pub voltage_v: f64,
    pub current_a: f64,
    pub energy_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn assemble(
    graph: &CircuitGraph,
    mutual: &[MutualInductance],
    nodes: &[usize],
    reference: &[usize],
    coordinates: &[usize],
    external_ports: usize,
    reference_ohm: f64,
) -> Result<CircuitDae> {
    let mut lowered = graph.clone();
    for (p, &terminal) in coordinates.iter().enumerate() {
        let a = nodes[terminal];
        let b = nodes[reference[terminal]];
        lowered.branches.push((
            a,
            b,
            Branch::Resistor {
                ohms: reference_ohm,
            },
        ));
        lowered.branches.push((
            a,
            b,
            Branch::CurrentSource {
                port: external_ports + p,
            },
        ));
    }
    Ok(assemble_circuit_with_mutual_inductance(&lowered, mutual)?)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    reference: &[usize],
    vertex_count: usize,
    graph: CircuitGraph,
    terminal_nodes: Vec<usize>,
    mutual: Vec<MutualInductance>,
    initial: SpatialTerminalCircuitInitial,
    config: SpatialTerminalCircuitConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Rc<Layout>, Vec<f64>)> {
    poll(cancelled)?;
    let count = reference.len();
    if terminal_nodes.len() != count
        || terminal_nodes.iter().any(|&n| n >= graph.node_count)
        || reference.iter().any(|&r| r >= count)
    {
        return Err(invalid("complete original electrode-to-node map required"));
    }
    let coordinates: Vec<_> = reference.iter().enumerate()
        // A real wire fixes zero relative voltage without a self-loop source.
        .filter_map(|(t, &r)| (terminal_nodes[t] != terminal_nodes[r]).then_some(t))
        .collect();
    let k = coordinates.len();
    if k > MAX_TERMINAL_CIRCUIT_COORDINATES
        || !(2..=MAX_TERMINAL_CIRCUIT_NODES).contains(&graph.node_count)
        || graph.branches.len().saturating_add(2 * k) > MAX_TERMINAL_CIRCUIT_BRANCHES
        || graph.transformers.len() > 16
        || mutual.len() > graph.branches.len().saturating_mul(graph.branches.len()) / 2
        || !positive(config.reference_resistance_ohm)
    {
        return Err(invalid(
            "bounded native graph, complete electrode nodes and positive reference required",
        ));
    }
    let mut source_kinds = vec![None; MAX_TERMINAL_CIRCUIT_SOURCE_PORTS];
    let (mut external_ports, mut inductors, mut capacitors, mut voltage_sources) = (0, 0, 0, 0);
    for &(_, _, branch) in &graph.branches {
        poll(cancelled)?;
        let (port, kind) = match branch {
            Branch::Resistor { ohms } => {
                if !positive(ohms) {
                    return Err(invalid("representable physical resistance required"));
                }
                continue;
            }
            Branch::Inductor { henries } => {
                if !positive(henries) {
                    return Err(invalid("representable physical inductance required"));
                }
                inductors += 1;
                continue;
            }
            Branch::Capacitor { farads } => {
                if !positive(farads) {
                    return Err(invalid("representable physical capacitance required"));
                }
                capacitors += 1;
                continue;
            }
            Branch::VoltageSource { port } => {
                voltage_sources += 1;
                (port, true)
            }
            Branch::CurrentSource { port } => (port, false),
        };
        if port >= source_kinds.len() || source_kinds[port].is_some_and(|old| old != kind) {
            return Err(invalid(
                "source port exceeds capacity or mixes voltage/current units",
            ));
        }
        source_kinds[port] = Some(kind);
        external_ports = external_ports.max(port + 1);
    }
    let n = inductors + 2 * capacitors + graph.node_count - 1
        + voltage_sources
        + graph.transformers.len();
    let payload = 16 * k * k + 12 * k + 8 * n + 6 * graph.branches.len() + 4 * external_ports;
    let payload = count
        .checked_mul(6)
        .and_then(|v| v.checked_add(payload))
        .and_then(|v| {
            vertex_count
                .checked_mul(3)
                .and_then(|vertices| v.checked_add(vertices))
        });
    if external_ports
        .checked_add(k)
        .is_none_or(|ports| ports > MAX_TERMINAL_CIRCUIT_PORTS)
        || n > MAX_TERMINAL_CIRCUIT_STATES
        || n > config.maximum_native_states
        || payload.is_none_or(|v| v > config.maximum_interface_scalars)
        || initial.fluxes_wb.len() != inductors
        || initial.charges_c.len() != capacitors
        || initial
            .fluxes_wb
            .iter()
            .chain(&initial.charges_c)
            .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "bounded native states and complete original circuit storage required",
        ));
    }
    poll(cancelled)?;
    let dae = assemble(
        &graph,
        &mutual,
        &terminal_nodes,
        reference,
        &coordinates,
        external_ports,
        config.reference_resistance_ohm,
    )?;
    poll(cancelled)?;
    if dae.ports != external_ports + k {
        return Err(invalid(
            "native assembly changed the original/private port map",
        ));
    }
    let layout = Rc::new(Layout {
        graph,
        mutual,
        dae,
        nodes: terminal_nodes,
        reference: reference.to_vec(),
        coordinates,
        source_kinds,
        external_ports,
        reference_ohm: config.reference_resistance_ohm,
    });
    admit_inputs(&layout, &initial.external_inputs)?;
    let mut state = vec![0.0; layout.dae.system.state_dim()];
    for (&index, &flux) in layout.dae.flux_index.iter().zip(&initial.fluxes_wb) {
        state[index] = flux;
    }
    for (&index, &charge) in layout.dae.charge_index.iter().zip(&initial.charges_c) {
        state[index] = charge;
    }
    let mut input = initial.external_inputs;
    input.resize(layout.dae.ports, 0.0);
    let state = native_seed(&layout.dae, &state, &input).map_err(|error| match error {
        NonlinearWindingError::Cancelled => SpatialTerminalCircuitError::Cancelled,
        other => SpatialTerminalCircuitError::Initial(other),
    })?;
    poll(cancelled)?;
    Ok((layout, state))
}

pub(super) fn advance(
    layout: &Layout,
    state: &[f64],
    input: &[f64],
    dt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StepRecord> {
    poll(cancelled)?;
    // A finite-stage source is not a new initial condition. Original accepted
    // flux/charge and the descriptor's own stage multipliers stay intact.
    let (step, _) = layout.dae.step_audited(state, input, dt)?;
    poll(cancelled)?;
    Ok(step)
}

pub(super) fn rebind_resistors(
    layout: &Layout,
    state: &[f64],
    sources: &[GraphResistorMaterial],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Rc<Layout>> {
    poll(cancelled)?;
    let count = layout
        .graph
        .branches
        .iter()
        .filter(|(_, _, branch)| matches!(branch, Branch::Resistor { .. }))
        .count();
    if sources.len() != count || sources.len() > MAX_TERMINAL_CIRCUIT_BRANCHES {
        return Err(invalid(
            "every physical graph resistor requires one bounded source",
        ));
    }
    let mut graph = layout.graph.clone();
    for (index, source) in sources.iter().enumerate() {
        poll(cancelled)?;
        if index > 0 && sources[index - 1].branch >= source.branch
            || !source
                .conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
            || !positive(source.conductor.resistance_ohm())
        {
            return Err(invalid(
                "unique original static Ohmic graph sources required",
            ));
        }
        let Some((_, _, branch @ Branch::Resistor { .. })) = graph.branches.get_mut(source.branch)
        else {
            return Err(invalid(
                "source does not name an original physical graph resistor",
            ));
        };
        *branch = Branch::Resistor {
            ohms: source.conductor.resistance_ohm(),
        };
    }
    let dae = assemble(
        &graph,
        &layout.mutual,
        &layout.nodes,
        &layout.reference,
        &layout.coordinates,
        layout.external_ports,
        layout.reference_ohm,
    )?;
    let old = &layout.dae;
    if dae.flux_index != old.flux_index
        || dae.charge_index != old.charge_index
        || dae.node_potential_index != old.node_potential_index
        || dae.ports != old.ports
        || dae.system.state_dim() != old.system.state_dim()
        || dae.system.hamiltonian(state).to_bits() != old.system.hamiltonian(state).to_bits()
    {
        return Err(invalid(
            "resistance refresh changed original native storage or terminal map",
        ));
    }
    poll(cancelled)?;
    Ok(Rc::new(Layout {
        graph,
        dae,
        mutual: layout.mutual.clone(),
        nodes: layout.nodes.clone(),
        reference: layout.reference.clone(),
        coordinates: layout.coordinates.clone(),
        source_kinds: layout.source_kinds.clone(),
        external_ports: layout.external_ports,
        reference_ohm: layout.reference_ohm,
    }))
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn account<F: FieldReceipt>(
    layout: &Layout,
    state: &[f64],
    native: StepRecord,
    field: F,
    graph_materials: &[GraphResistorMaterial],
    input: &[f64],
    external: &[f64],
    dt: f64,
    controls: Allowances,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialTerminalCircuitStep<F>> {
    let account = field.account();
    let dae = &layout.dae;
    if account.terminals.len() != layout.nodes.len()
        || account.terminals.iter().enumerate().any(|(t, transfer)| {
            transfer.terminal != t || transfer.reference_terminal != layout.reference[t]
        })
    {
        return Err(invalid(
            "field receipt lost its original physical electrode order",
        ));
    }
    let midpoint = state
        .iter()
        .zip(&native.x)
        .map(|(&a, &b)| f64::midpoint(a, b))
        .collect::<Vec<_>>();
    let mut nodes = vec![0.0; layout.graph.node_count];
    for (node, &index) in dae.node_potential_index.iter().enumerate() {
        nodes[node + 1] = finite(midpoint[index])?;
    }
    let (mut maximum_voltage, mut maximum_current, mut reference_heat, mut interface_work) =
        (0.0_f64, 0.0_f64, 0.0, 0.0);
    for (p, &terminal) in layout.coordinates.iter().enumerate() {
        poll(cancelled)?;
        let reference = layout.reference[terminal];
        let v = finite(nodes[layout.nodes[terminal]] - nodes[layout.nodes[reference]])?;
        let current = finite(v / layout.reference_ohm - input[layout.external_ports + p])?;
        maximum_voltage = maximum_voltage
            .max(finite(v - account.terminals[terminal].voltage_to_reference_v)?.abs());
        maximum_current = maximum_current
            .max(finite(current - account.terminals[terminal].incoming_current_a)?.abs());
        close(
            v - native.y[layout.external_ports + p],
            controls.voltage_v,
            "native private port and original node voltage disagree",
        )?;
        reference_heat = finite(reference_heat + dt * v * (v / layout.reference_ohm))?;
        interface_work = finite(
            interface_work
                + dt * input[layout.external_ports + p] * native.y[layout.external_ports + p],
        )?;
    }
    close(
        maximum_voltage,
        controls.voltage_v,
        "native circuit and field electrode voltages disagree",
    )?;
    close(
        maximum_current,
        controls.current_a,
        "native circuit and field electrode currents disagree",
    )?;
    close(
        interface_work + account.terminal_work_j - reference_heat,
        controls.energy_j,
        "private circuit work does not cancel physical terminal work and reference loss",
    )?;

    let currents = dae.inductor_currents(&midpoint)?;
    let (mut l, mut c, mut source) = (0, 0, 0);
    let source_start =
        dae.system.differential_dim() + layout.graph.node_count - 1 + dae.charge_index.len();
    let mut kcl = vec![0.0; layout.graph.node_count];
    let mut branches = Vec::with_capacity(layout.graph.branches.len());
    let mut graph_resistors = Vec::new();
    let (mut resistor_heat, mut circuit_change) = (0.0, 0.0);
    for (branch, &(a, b, kind)) in layout.graph.branches.iter().enumerate() {
        poll(cancelled)?;
        let v = finite(nodes[a] - nodes[b])?;
        let current = match kind {
            Branch::Resistor { ohms } => finite(v / ohms)?,
            Branch::Inductor { .. } => {
                let index = dae.flux_index[l];
                let current = currents[l];
                circuit_change =
                    finite(circuit_change + current * finite(native.x[index] - state[index])?)?;
                close(
                    finite(native.x[index] - state[index])? - dt * v,
                    dt * controls.voltage_v,
                    "original inductor flux increment violates terminal voltage",
                )?;
                l += 1;
                current
            }
            Branch::Capacitor { farads } => {
                let index = dae.charge_index[c];
                let dq = finite(native.x[index] - state[index])?;
                let effort = finite(midpoint[index] / farads)?;
                close(
                    v - effort,
                    controls.voltage_v,
                    "original capacitor charge violates terminal voltage",
                )?;
                circuit_change = finite(circuit_change + effort * dq)?;
                c += 1;
                finite(dq / dt)?
            }
            Branch::VoltageSource { port } => {
                close(
                    v - external[port],
                    controls.voltage_v,
                    "original voltage source does not hold its authored effort",
                )?;
                let value = finite(midpoint[source_start + source])?;
                source += 1;
                value
            }
            Branch::CurrentSource { port } => -external[port],
        };
        let power = finite(v * current)?;
        let work = finite(dt * power)?;
        if matches!(kind, Branch::Resistor { .. }) {
            if power < 0.0 || work < 0.0 || current != 0.0 && (power == 0.0 || work == 0.0) {
                return Err(invalid("physical graph resistor heat is not representable"));
            }
            resistor_heat = finite(resistor_heat + work)?;
            graph_resistors.push(WindingGraphResistorHeat {
                branch,
                voltage_v: v,
                current_a: current,
                power_w: power,
                heat_j: work,
            });
        }
        kcl[a] = finite(kcl[a] + current)?;
        kcl[b] = finite(kcl[b] - current)?;
        branches.push(SpatialTerminalCircuitBranch {
            branch,
            voltage_v: v,
            current_a: current,
            absorbed_work_j: work,
        });
    }
    let transformer_start = dae.system.state_dim() - layout.graph.transformers.len();
    for (index, transformer) in layout.graph.transformers.iter().enumerate() {
        poll(cancelled)?;
        let primary = finite(midpoint[transformer_start + index])?;
        let secondary = finite(-transformer.ratio * primary)?;
        let (a, b) = transformer.primary;
        let (c, d) = transformer.secondary;
        let vp = finite(nodes[a] - nodes[b])?;
        let vs = finite(nodes[c] - nodes[d])?;
        close(
            vp - transformer.ratio * vs,
            controls.voltage_v,
            "original ideal transformer voltage constraint failed",
        )?;
        close(
            dt * (vp * primary + vs * secondary),
            controls.energy_j,
            "original ideal transformer creates work",
        )?;
        for ((a, b), current) in [
            (transformer.primary, primary),
            (transformer.secondary, secondary),
        ] {
            kcl[a] = finite(kcl[a] + current)?;
            kcl[b] = finite(kcl[b] - current)?;
        }
    }
    let mut component_balance = std::collections::BTreeMap::new();
    for (terminal, transfer) in account.terminals.iter().enumerate() {
        poll(cancelled)?;
        let node = layout.nodes[terminal];
        kcl[node] = finite(kcl[node] + transfer.incoming_current_a)?;
        let balance = component_balance.entry(transfer.component).or_insert(0.0);
        *balance = finite(*balance + transfer.incoming_current_a)?;
        // Also check genuinely shorted electrodes and the reference itself.
        let actual = finite(nodes[node] - nodes[layout.nodes[layout.reference[terminal]]])?;
        close(
            actual - transfer.voltage_to_reference_v,
            controls.voltage_v,
            "physical electrode voltage disagrees with its original graph nodes",
        )?;
    }
    for balance in component_balance.into_values() {
        close(
            balance,
            controls.current_a,
            "physical conductor electrode currents do not balance",
        )?;
    }
    let maximum_kcl = kcl.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    close(
        maximum_kcl,
        controls.current_a,
        "original graph and physical electrodes violate Kirchhoff current balance",
    )?;
    let mut external_work = 0.0;
    for (&u, &y) in external.iter().zip(&native.y) {
        external_work = finite(external_work + dt * u * y)?;
    }
    close(
        native.delta_h - circuit_change,
        controls.energy_j,
        "native and original constitutive storage changes disagree",
    )?;
    close(
        native.dissipated - resistor_heat - reference_heat,
        controls.energy_j,
        "physical resistor and canceled reference heat do not reproduce native dissipation",
    )?;
    close(
        native.supplied - external_work - interface_work,
        controls.energy_j,
        "original external source and private interface work disagree with native supply",
    )?;
    close(
        native.balance_residual(),
        controls.energy_j,
        "original native circuit energy does not close",
    )?;
    let circuit_before = finite(dae.system.hamiltonian(state))?;
    let circuit_after = finite(dae.system.hamiltonian(&native.x))?;
    let stored_before_j = finite(account.stored_before_j + circuit_before)?;
    let stored_after_j = finite(account.stored_after_j + circuit_after)?;
    let stored_change_j = finite(account.stored_change_j + circuit_change)?;
    let heat_j = finite(account.heat_j + resistor_heat)?;
    let numerical_dissipation_j = finite(account.numerical_dissipation_j)?;
    let impressed_work_j = finite(account.impressed_work_j)?;
    let supplied_work_j = finite(external_work + impressed_work_j)?;
    let parameter_work_j = finite(account.parameter_work_j)?;
    let mechanical_work_j = finite(account.mechanical_work_j)?;
    if circuit_before < 0.0 || circuit_after < 0.0 || heat_j < 0.0 || numerical_dissipation_j < 0.0
    {
        return Err(invalid(
            "physical circuit/field storage and losses violate passivity",
        ));
    }
    let energy_residual = finite(
        stored_change_j + heat_j + numerical_dissipation_j + mechanical_work_j
            - supplied_work_j
            - parameter_work_j,
    )?;
    close(
        energy_residual,
        controls.energy_j,
        "coupled circuit/field storage, physical heat, numerical loss and source work do not close",
    )?;
    let terminal_voltage_v = layout.nodes.iter().map(|&n| nodes[n]).collect();
    let terminal_current_a = account
        .terminals
        .iter()
        .map(|t| t.incoming_current_a)
        .collect();
    let external_source_outputs = native.y[..layout.external_ports].to_vec();
    poll(cancelled)?;
    Ok(SpatialTerminalCircuitStep {
        native,
        field,
        node_voltages_v: nodes,
        terminal_voltage_v,
        terminal_current_a,
        branches,
        graph_resistors,
        graph_materials: graph_materials.to_vec(),
        external_source_outputs,
        resistor_heat_j: resistor_heat,
        cancelled_reference_heat_j: reference_heat,
        circuit_stored_before_j: circuit_before,
        circuit_stored_after_j: circuit_after,
        circuit_stored_change_j: circuit_change,
        stored_before_j,
        stored_after_j,
        stored_change_j,
        heat_j,
        numerical_dissipation_j,
        supplied_work_j,
        impressed_work_j,
        parameter_work_j,
        mechanical_work_j,
        external_work_j: external_work,
        energy_residual_j: energy_residual,
        relative_terminal_residual: 0.0,
        maximum_terminal_voltage_defect_v: maximum_voltage,
        maximum_terminal_current_defect_a: maximum_current,
        maximum_kirchhoff_current_defect_a: maximum_kcl,
        field_trials: 0,
        reused_response: false,
        report: None,
    })
}
