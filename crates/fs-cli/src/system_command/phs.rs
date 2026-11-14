//! A linear port-Hamiltonian system written as matrices: the exact
//! held-input owner of `fs-system`, with its heat outlet.

use core::num::NonZeroUsize;

use fs_ports::StableId;
use fs_ports::component::{PortRole, PortSpec};
use fs_system::phs_exact::{
    ExactPortHamiltonianComponent, LinearSplitBudget, PortHamiltonian, QuadraticStorage,
};

use super::{Domain, J, Refusal, coefficient_matrix, fields};

/// State coordinates one element may hold: the owner's dense preparation
/// capacity, which counts every coordinate.
pub(super) const MAX_STATES: usize = 64;
/// Power ports one element may have.
const MAX_PORTS: usize = 64;

/// A `port_hamiltonian`: `dx/dt = (J - R) Q x + G u` with stored energy
/// `H = x' Q x / 2` and port outputs `y = G' Q x`. Each port receives one
/// input `u` and returns the mean of its `y` over the step; the owner steps
/// the held input exactly and sends the dissipated work out through `heat`.
/// Returns the owner and its initial state.
pub(super) fn port_hamiltonian(
    definition: &J,
    what: &str,
) -> Result<(ExactPortHamiltonianComponent, Vec<f64>), Refusal> {
    fields(
        definition,
        &[
            "kind",
            "energy",
            "structure",
            "dissipation",
            "input",
            "ports",
            "state",
        ],
        what,
    )?;
    let states = definition
        .get("energy")
        .and_then(J::as_array)
        .map(<[J]>::len)
        .filter(|rows| (1..=MAX_STATES).contains(rows))
        .ok_or_else(|| format!("{what}.energy: Q must have 1 to {MAX_STATES} square rows"))?;
    let declared = definition
        .get("ports")
        .and_then(J::as_array)
        .filter(|ports| (1..=MAX_PORTS).contains(&ports.len()))
        .ok_or_else(|| format!("{what}.ports: name 1 to {MAX_PORTS} power ports"))?;
    let ports = declared
        .iter()
        .enumerate()
        .map(|(index, port)| power_port(port, &format!("{what}.ports[{index}]")))
        .collect::<Result<Vec<_>, _>>()?;
    let matrix = |key: &str, columns: usize| {
        coefficient_matrix(
            definition.get(key),
            states,
            columns,
            &format!("{what}.{key}"),
        )
    };
    let energy = matrix("energy", states)?;
    let structure = matrix("structure", states)?;
    let dissipation = matrix("dissipation", states)?;
    let input = matrix("input", ports.len())?;
    let state = match definition.get("state") {
        None => vec![0.0; states],
        Some(value) => value
            .as_array()
            .filter(|values| values.len() == states)
            .and_then(|values| {
                values
                    .iter()
                    .map(|value| value.as_f64().filter(|value| value.is_finite()))
                    .collect::<Option<Vec<_>>>()
            })
            .ok_or_else(|| format!("{what}.state: {states} finite numbers, one per row of Q"))?,
    };
    let storage =
        QuadraticStorage::new(energy, states).map_err(|e| format!("{what}.energy: {e}"))?;
    let system = PortHamiltonian::new(
        states,
        ports.len(),
        structure,
        dissipation,
        input,
        Box::new(storage),
    )
    .map_err(|e| format!("{what}: {e}"))?;
    let budget = LinearSplitBudget {
        max_modes: MAX_STATES,
        max_omega_dt: 1.0,
    };
    let owner = ExactPortHamiltonianComponent::new(system, ports, budget)
        .map_err(|e| format!("{what}: {e}"))?;
    Ok((owner, state))
}

fn power_port(port: &J, what: &str) -> Result<PortSpec, String> {
    fields(port, &["name", "domain", "role"], what)?;
    let name = port
        .str_field("name")
        .ok_or_else(|| format!("{what}: missing name"))?;
    let id = StableId::new(name)
        .map_err(|_| format!("{what}: name {name:?} is not a port identifier"))?;
    let domain = port
        .str_field("domain")
        .and_then(Domain::parse)
        .ok_or_else(|| format!("{what}: domain is one of translation, rotation, electrical"))?;
    let role = match port.str_field("role") {
        Some("effort_in") => PortRole::EffortIn,
        Some("flow_in") => PortRole::FlowIn,
        _ => {
            return Err(format!(
                "{what}: role is \"effort_in\" (receives an effort, returns its mean \
                 flow) or \"flow_in\" (receives a flow, returns its mean effort)"
            ));
        }
    };
    PortSpec::power(id, domain.port_kind(), role, NonZeroUsize::MIN)
        .map_err(|e| format!("{what}: {e:?}"))
}
