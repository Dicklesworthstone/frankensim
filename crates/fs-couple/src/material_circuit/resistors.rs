//! Original source lowering and midpoint branch heat, shared by circuit consumers.
//! These helpers neither solve electrical storage nor infer heat from a ledger.

use super::*;

/// Replace every resistor slot with its complete original conductor source.
/// All circuit consumers use the same branch-count, geometry and scale gates.
pub(crate) fn lower_resistors(
    graph: &CircuitGraph,
    mut sources: Vec<MaterialResistor>,
) -> Result<(CircuitGraph, Vec<MaterialResistor>), MaterialCircuitError> {
    use MaterialCircuitError::Invalid;
    let count = admit_graph(graph)?;
    if sources.len() != count {
        return Err(Invalid(
            "every resistor requires exactly one material binding",
        ));
    }
    sources.sort_by_key(|source| source.branch);
    let mut lowered = graph.clone();
    for (index, source) in sources.iter().enumerate() {
        if index > 0 && sources[index - 1].branch == source.branch {
            return Err(Invalid("duplicate resistor binding"));
        }
        let branch = lowered
            .branches
            .get_mut(source.branch)
            .ok_or(Invalid("material branch index out of range"))?;
        if !matches!(branch.2, Branch::Resistor { .. }) {
            return Err(Invalid("material binding does not name a resistor"));
        }
        let resistance = source.conductor.resistance_ohm();
        if !positive(resistance) || !positive(1.0 / resistance) {
            return Err(Invalid("resistor conductance is unrepresentable"));
        }
        branch.2 = Branch::Resistor { ohms: resistance };
    }
    Ok((lowered, sources))
}

/// Independently evaluate original material/contact power at native midpoint
/// node efforts. Whole-circuit storage and other loss channels are checked by
/// their owner after adding this explicit branch-local heat exactly once.
#[allow(clippy::too_many_arguments)]
pub(crate) fn account_resistor_heat(
    graph: &CircuitGraph,
    resistors: &[MaterialResistor],
    contacts: &[SeriesElectricalContact],
    node_indices: &[usize],
    x0: &[f64],
    record: &StepRecord,
    dt: f64,
    controls: MaterialCircuitSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<ResistorHeat>, Vec<ContactHeat>, f64), MaterialCircuitError> {
    let voltage = |node: usize| {
        if node == 0 {
            0.0
        } else {
            let index = node_indices[node - 1];
            f64::midpoint(x0[index], record.x[index])
        }
    };
    let mut heat_j = 0.0;
    let mut heat = Vec::with_capacity(resistors.len());
    let mut contact_heat = Vec::with_capacity(contacts.len());
    for source in resistors {
        poll(cancelled)?;
        let (a, b, branch) = graph.branches[source.branch];
        let terminal_v = voltage(a) - voltage(b);
        let has_contact = contacts
            .iter()
            .any(|contact| contact.branch() == source.branch);
        // Preserve the original contact-free conductor arithmetic.
        let (v, i) = if has_contact {
            let Branch::Resistor { ohms } = branch else {
                unreachable!("admitted resistor")
            };
            let current = terminal_v / ohms;
            if !current.is_finite() || (terminal_v != 0.0 && current == 0.0) {
                return Err(MaterialCircuitError::Balance(
                    "series current is unrepresentable",
                ));
            }
            (current * source.conductor.resistance_ohm(), current)
        } else {
            (
                terminal_v,
                source
                    .conductor
                    .current_a_for_drive(OhmicDrive::Voltage(terminal_v))?,
            )
        };
        let power = source.conductor.joule_power_w(i)?;
        let joules = dt * power;
        let defect = v * i - power;
        if !joules.is_finite()
            || (power != 0.0 && joules == 0.0)
            || !defect.is_finite()
            || defect.abs() > controls.power_tolerance_w
        {
            return Err(MaterialCircuitError::Balance(
                "resistor electrical work and material heat disagree",
            ));
        }
        heat_j += joules;
        let mut branch_power = power;
        for contact in contacts
            .iter()
            .filter(|contact| contact.branch() == source.branch)
        {
            poll(cancelled)?;
            let loss = contact.heat(i, dt)?;
            branch_power += loss.power_w;
            heat_j += loss.heat_j;
            contact_heat.push(loss);
        }
        let terminal_defect = terminal_v * i - branch_power;
        if !terminal_defect.is_finite() || terminal_defect.abs() > controls.power_tolerance_w {
            return Err(MaterialCircuitError::Balance(
                "bulk and contact powers do not close terminal work",
            ));
        }
        heat.push(ResistorHeat {
            material: source.clone(),
            voltage_v: v,
            current_a: i,
            power_w: power,
            heat_j: joules,
        });
    }
    contact_heat.sort_by_key(|heat| heat.contact.id());
    if !heat_j.is_finite() {
        return Err(MaterialCircuitError::Balance(
            "summed resistor heat is not representable",
        ));
    }
    Ok((heat, contact_heat, heat_j))
}
