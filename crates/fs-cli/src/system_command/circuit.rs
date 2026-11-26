//! `kirchhoff_circuit`: fs-system's [`KirchhoffCircuit`], a netlist of
//! `fs_phs::circuit` (resistors, inductors with mutual coupling, capacitors,
//! sources driven by external ports and ideal transformers) as one element
//! with a port per source and a heat port per resistor.
//!
//! The scene authors the graph, the start and the energy tolerance; every
//! step is the owner's own. The resistances are constant: the Component's
//! resistance laws and its operating-point iteration are not exposed.

use fs_ports::component::{PortRole, PortSpec};
use fs_system::circuit::{
    Branch, CircuitGraph, KirchhoffCircuit, MutualInductance, TransformerLink,
};
use fs_system::{NodeId, System};

use super::{J, Refusal, fields, index_below, number, number_or, positive};

/// The most external ports a scene's circuit may name: 0 to 63.
const PORTS: usize = 64;
/// The most nodes a scene's circuit may have, ground included.
const NODES: usize = 4096;
/// The field that makes a branch what it is.
const ELEMENTS: [&str; 5] = [
    "resistance_ohm",
    "inductance_h",
    "capacitance_f",
    "voltage_port",
    "current_port",
];

/// What a report reads besides the ports: the name of each state coordinate
/// in the owner's order (every inductor's flux, then every capacitor's
/// charge, each in branch order), which also names the adaptive accounts,
/// and the tolerance the run used.
#[derive(Clone, Debug)]
pub(super) struct Declared {
    pub(super) coordinates: Vec<String>,
    energy_tolerance_j: f64,
}

/// A whole number below `limit` at `key`.
fn whole(value: &J, key: &str, limit: usize, what: &str) -> Result<usize, String> {
    value
        .get(key)
        .and_then(|value| index_below(value, limit))
        .ok_or_else(|| format!("{what}.{key}: a whole number below {limit}"))
}

/// Two whole numbers below `limit` at `key`.
fn pair(value: &J, key: &str, limit: usize, what: &str) -> Result<(usize, usize), String> {
    match value.get(key).and_then(J::as_array) {
        Some([a, b]) => index_below(a, limit).zip(index_below(b, limit)),
        _ => None,
    }
    .ok_or_else(|| format!("{what}.{key}: two whole numbers below {limit}"))
}

/// The list at `key`, empty when the scene gives none.
fn list<'j>(value: &'j J, key: &str, what: &str) -> Result<&'j [J], String> {
    match value.get(key) {
        None => Ok(&[]),
        Some(list) => list
            .as_array()
            .ok_or_else(|| format!("{what}.{key}: a list")),
    }
}

/// A `kirchhoff_circuit` and its start: each inductor's `flux_wb` and each
/// capacitor's `charge_c` (0 when not given), in the owner's order.
pub(super) fn circuit(
    definition: &J,
    what: &str,
) -> Result<(KirchhoffCircuit, Vec<f64>, Declared), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "node_count",
            "branches",
            "mutual",
            "transformers",
            "energy_tolerance_j",
        ],
        what,
    )?;
    let node_count = whole(definition, "node_count", NODES + 1, what)?;
    let rows = definition
        .get("branches")
        .and_then(J::as_array)
        .ok_or_else(|| format!("{what}.branches: a list of branches"))?;
    let mut branches = Vec::with_capacity(rows.len());
    let (mut fluxes, mut charges) = (Vec::new(), Vec::new());
    let (mut flux_names, mut charge_names) = (Vec::new(), Vec::new());
    for (index, row) in rows.iter().enumerate() {
        let what = format!("{what}.branches[{index}]");
        let mut allowed = vec!["nodes", "flux_wb", "charge_c"];
        allowed.extend(ELEMENTS);
        fields(row, &allowed, &what)?;
        let (a, b) = pair(row, "nodes", node_count, &what)?;
        let given: Vec<&str> = ELEMENTS
            .into_iter()
            .filter(|key| row.get(key).is_some())
            .collect();
        let &[key] = given.as_slice() else {
            return Err(format!("{what}: exactly one of {}", ELEMENTS.join(", ")).into());
        };
        let element = match key {
            "resistance_ohm" => Branch::Resistor {
                ohms: number(row, key, &what)?,
            },
            "inductance_h" => {
                fluxes.push(number_or(row, "flux_wb", 0.0, &what)?);
                flux_names.push(format!("flux{index}_wb"));
                Branch::Inductor {
                    henries: number(row, key, &what)?,
                }
            }
            "capacitance_f" => {
                charges.push(number_or(row, "charge_c", 0.0, &what)?);
                charge_names.push(format!("charge{index}_c"));
                Branch::Capacitor {
                    farads: number(row, key, &what)?,
                }
            }
            "voltage_port" => Branch::VoltageSource {
                port: whole(row, key, PORTS, &what)?,
            },
            _ => Branch::CurrentSource {
                port: whole(row, key, PORTS, &what)?,
            },
        };
        for (start, owner) in [("flux_wb", "inductance_h"), ("charge_c", "capacitance_f")] {
            if row.get(start).is_some() && key != owner {
                return Err(format!("{what}: {start} belongs to a branch with {owner}").into());
            }
        }
        branches.push((a, b, element));
    }
    let (mutual, transformers) = couplings(definition, branches.len(), node_count, what)?;
    let energy_tolerance_j = positive(definition, "energy_tolerance_j", what)?;
    let graph = CircuitGraph {
        node_count,
        branches,
        transformers,
    };
    let circuit = KirchhoffCircuit::new(graph, mutual, energy_tolerance_j)
        .map_err(|e| format!("{what}: {e}"))?;
    let state = circuit
        .state(&fluxes, &charges)
        .map_err(|e| format!("{what}: {e}"))?;
    flux_names.extend(charge_names);
    let declared = Declared {
        coordinates: flux_names,
        energy_tolerance_j,
    };
    Ok((circuit, state, declared))
}

/// The mutual inductances between the `branch_count` branches and the ideal
/// transformers between the `node_count` nodes.
fn couplings(
    definition: &J,
    branch_count: usize,
    node_count: usize,
    what: &str,
) -> Result<(Vec<MutualInductance>, Vec<TransformerLink>), String> {
    let mutual = list(definition, "mutual", what)?
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let what = format!("{what}.mutual[{index}]");
            fields(row, &["branches", "inductance_h"], &what)?;
            Ok(MutualInductance {
                branches: pair(row, "branches", branch_count, &what)?,
                henries: number(row, "inductance_h", &what)?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let transformers = list(definition, "transformers", what)?
        .iter()
        .enumerate()
        .map(|(index, row)| {
            let what = format!("{what}.transformers[{index}]");
            fields(row, &["primary", "secondary", "ratio"], &what)?;
            Ok(TransformerLink {
                primary: pair(row, "primary", node_count, &what)?,
                secondary: pair(row, "secondary", node_count, &what)?,
                ratio: number(row, "ratio", &what)?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok((mutual, transformers))
}

/// The state, one value per coordinate name; then each source port's last
/// output (the current a voltage port delivers into the circuit, or the
/// voltage a current port drives at), each resistor's last heat rate by its
/// branch index (null before the first step), and the tolerance the run
/// used.
pub(super) fn values(
    system: &System,
    node: NodeId,
    ports: &[PortSpec],
    declared: &Declared,
) -> Vec<(String, f64)> {
    let state = system.state::<Vec<f64>>(node);
    let output = |port: &str| {
        system
            .port(node, port)
            .filter(|_| system.steps() > 0)
            .and_then(|port| system.last_output(port).first().copied())
            .unwrap_or(f64::NAN)
    };
    let mut values: Vec<(String, f64)> = declared
        .coordinates
        .iter()
        .enumerate()
        .map(|(index, name)| {
            let value = state.and_then(|state| state.get(index)).copied();
            (name.clone(), value.unwrap_or(f64::NAN))
        })
        .collect();
    for port in ports {
        let id = port.id().as_str();
        if let Some(branch) = id.strip_prefix("heat") {
            values.push((format!("heat{branch}_w"), -output(id)));
        } else {
            let returned = if port.role() == PortRole::EffortIn {
                "current_a"
            } else {
                "voltage_v"
            };
            values.push((format!("{id}.{returned}"), output(id)));
        }
    }
    values.push(("energy_tolerance_j".to_owned(), declared.energy_tolerance_j));
    values
}
