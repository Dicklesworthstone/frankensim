//! Close the actual source and external loads before requiring an open-circuit
//! impedance. Component or open-port poles need not be poles of the loaded graph.

use super::*;

/// Solve simultaneous phase-bearing volume-flow sources and explicit passive
/// loads at one frequency. The sources are complex peak amplitudes, all on the
/// same exp(-i omega t) clock. Every non-driven port belongs to exactly one
/// scalar load or the optional coupled boundary. A zero source remains an
/// imposed zero-flow port; no source is replaced by a passive termination.
///
/// One factorization serves the complete source vector, including at an
/// open-circuit pole regularized by the actual loads. All reported work uses
/// the combined fields and therefore includes the interference cross terms.
pub fn network_drive_coherent(
    net: &DuctNetwork,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    sources: &[PortFlowSource],
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkCoherentResponse, DuctError> {
    Ok(network_drive_sources(
        net,
        gas,
        omega,
        loss,
        &flow_sources(sources)?,
        loads,
        coupled,
    )?
    .fields)
}

/// [`network_drive_coherent`] with the original regional edge gases and
/// independent exterior. Relative source phases are not reset at interfaces.
pub fn network_drive_coherent_regional(
    net: &DuctNetwork,
    media: &NetworkMedia,
    omega: f64,
    loss: LossModel,
    sources: &[PortFlowSource],
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkCoherentResponse, DuctError> {
    Ok(network_drive_sources_regional(
        net,
        media,
        omega,
        loss,
        &flow_sources(sources)?,
        loads,
        coupled,
    )?
    .fields)
}

/// Solve coherent pressure generators with finite source impedance and/or
/// ideal prescribed flow sources. Pressure boundaries use `p + Zs U = p_open`;
/// every actual source flow therefore responds to the other sources, the graph
/// and the full exterior load. All phasors are peak amplitudes at this frequency.
/// A zero pressure drive retains its source impedance as a physical boundary.
///
/// Sources and passive loads partition every external port exactly once. Source
/// impedances must be finite with Re Zs >= 0; zero admits ideal pressure drive.
/// The loaded scattering solve remains available at an open-circuit pole.
pub fn network_drive_sources(
    net: &DuctNetwork,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    sources: &[PortSource],
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkSourceResponse, DuctError> {
    drive_sources_bound(net, gas, None, omega, loss, sources, loads, coupled)
}

/// [`network_drive_sources`] retaining the original per-edge and exterior gas
/// records. Every source pressure/flow phase refers to the same physical clock.
pub fn network_drive_sources_regional(
    net: &DuctNetwork,
    media: &NetworkMedia,
    omega: f64,
    loss: LossModel,
    sources: &[PortSource],
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkSourceResponse, DuctError> {
    media.validate_for(net)?;
    drive_sources_bound(
        net,
        media.exterior_gas(),
        Some(media),
        omega,
        loss,
        sources,
        loads,
        coupled,
    )
}

#[allow(clippy::too_many_arguments)]
fn drive_sources_bound(
    net: &DuctNetwork,
    exterior: &GasState,
    media: Option<&NetworkMedia>,
    omega: f64,
    loss: LossModel,
    sources: &[PortSource],
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkSourceResponse, DuctError> {
    let ordinary = network_port_impedance_bound(net, exterior, media, omega, loss)
        .and_then(|response| response.drive_sources(sources, loads, coupled));
    match ordinary {
        Ok(response) => Ok(response),
        Err(DuctError::Singular) => {
            let fields =
                drive_retained_sources(net, exterior, media, omega, loss, sources, loads, coupled)?;
            let labels: Vec<&str> = fields.labels.iter().map(String::as_str).collect();
            NetworkSourceResponse::from_fields(
                omega,
                &labels,
                &fields.conditions,
                fields.pressures,
                fields.inflows,
            )
        }
        Err(error) => Err(error),
    }
}

/// Unit volume flow into one active port with every other port explicitly
/// terminated. Unlike an open-circuit N-port query followed by `drive`, this
/// remains available when the supplied loads regularize an open-circuit pole.
/// The returned pressure/current order follows the original graph's port order.
pub fn network_drive(
    net: &DuctNetwork,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    active: &str,
    loads: &[PortTermination],
) -> Result<NetworkDriveResponse, DuctError> {
    drive_bound(net, gas, None, omega, loss, active, loads, None)
}

/// [`network_drive`] with an explicit passive mutual-impedance boundary.
/// Every non-active port belongs to exactly one scalar or coupled load.
pub fn network_drive_coupled(
    net: &DuctNetwork,
    gas: &GasState,
    omega: f64,
    loss: LossModel,
    active: &str,
    loads: &[PortTermination],
    coupled: &CoupledTermination,
) -> Result<NetworkDriveResponse, DuctError> {
    drive_bound(net, gas, None, omega, loss, active, loads, Some(coupled))
}

/// [`network_drive`] with the original edge gas records and exterior medium.
/// All radiation, chimney and porous-load admission uses those same records.
pub fn network_drive_regional(
    net: &DuctNetwork,
    media: &NetworkMedia,
    omega: f64,
    loss: LossModel,
    active: &str,
    loads: &[PortTermination],
) -> Result<NetworkDriveResponse, DuctError> {
    media.validate_for(net)?;
    drive_bound(
        net,
        media.exterior_gas(),
        Some(media),
        omega,
        loss,
        active,
        loads,
        None,
    )
}

/// [`network_drive_regional`] closed by a full passive exterior load matrix.
/// A singular load matrix is admissible: its constitutive equations are
/// retained, so pressure-release combinations need no inverse impedance.
pub fn network_drive_coupled_regional(
    net: &DuctNetwork,
    media: &NetworkMedia,
    omega: f64,
    loss: LossModel,
    active: &str,
    loads: &[PortTermination],
    coupled: &CoupledTermination,
) -> Result<NetworkDriveResponse, DuctError> {
    media.validate_for(net)?;
    drive_bound(
        net,
        media.exterior_gas(),
        Some(media),
        omega,
        loss,
        active,
        loads,
        Some(coupled),
    )
}

#[allow(clippy::too_many_arguments)]
fn drive_bound(
    net: &DuctNetwork,
    exterior: &GasState,
    media: Option<&NetworkMedia>,
    omega: f64,
    loss: LossModel,
    active: &str,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkDriveResponse, DuctError> {
    // Preserve the existing, smaller solve and its exact operation order
    // wherever the requested open-circuit elimination is admissible.
    let ordinary = network_port_impedance_bound(net, exterior, media, omega, loss)
        .and_then(|response| response.drive_with(active, loads, coupled, exterior));
    match ordinary {
        Ok(response) => Ok(response),
        Err(DuctError::Singular) => {
            drive_retained(net, exterior, media, omega, loss, active, loads, coupled)
        }
        Err(error) => Err(error),
    }
}

#[allow(clippy::too_many_arguments)]
fn drive_retained(
    net: &DuctNetwork,
    exterior: &GasState,
    media: Option<&NetworkMedia>,
    omega: f64,
    loss: LossModel,
    active_label: &str,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkDriveResponse, DuctError> {
    let sources = [PortSource {
        label: active_label.to_owned(),
        excitation: PortExcitation::VolumeFlow(C64::ONE),
    }];
    let fields =
        drive_retained_sources(net, exterior, media, omega, loss, &sources, loads, coupled)?;
    let active = fields
        .labels
        .iter()
        .position(|label| label == active_label)
        .expect("source port was admitted");
    Ok(NetworkDriveResponse {
        impedance: fields.pressures[active],
        port_pressures: fields.pressures,
        port_inflows: fields.inflows,
    })
}

struct LoadedFields {
    labels: Vec<String>,
    conditions: Vec<Option<PortExcitation>>,
    pressures: Vec<C64>,
    inflows: Vec<C64>,
}

#[allow(clippy::too_many_arguments)]
fn drive_retained_sources(
    net: &DuctNetwork,
    exterior: &GasState,
    media: Option<&NetworkMedia>,
    omega: f64,
    loss: LossModel,
    sources: &[PortSource],
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<LoadedFields, DuctError> {
    net.validate()?;
    validate_medium(exterior, omega)?;
    let port_nodes: Vec<usize> = net
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(node, value)| {
            matches!(value.boundary, NodeBoundary::Port { .. }).then_some(node)
        })
        .collect();
    let labels: Vec<&str> = port_nodes
        .iter()
        .map(|&node| match &net.nodes[node].boundary {
            NodeBoundary::Port { label } => label.as_str(),
            _ => unreachable!("port indices were selected above"),
        })
        .collect();
    let conditions = source_conditions(&labels, sources)?;
    if loads.len() + coupled.map_or(0, |value| value.labels.len()) != labels.len() - sources.len() {
        return Err(bad(
            "every other network port needs an explicit termination",
        ));
    }
    let mut resolved: Vec<Option<C64>> = conditions
        .iter()
        .map(|condition| condition.and_then(PortExcitation::pressure).map(|(_, z)| z))
        .collect();
    let mut coupled_rows = vec![None; labels.len()];
    let mut seen: Vec<bool> = conditions.iter().map(Option::is_some).collect();
    if let Some(value) = coupled {
        admit_coupled(value, omega)?;
        for (row, label) in value.labels.iter().enumerate() {
            let port = labels
                .iter()
                .position(|candidate| *candidate == label.as_str())
                .ok_or_else(|| bad("unknown network coupled-termination label"))?;
            if seen[port] {
                return Err(bad("duplicate or active-port coupled termination"));
            }
            seen[port] = true;
            coupled_rows[port] = Some(row);
        }
    }
    for load in loads {
        let port = labels
            .iter()
            .position(|candidate| *candidate == load.label.as_str())
            .ok_or_else(|| bad("unknown network termination label"))?;
        if seen[port] {
            return Err(bad("duplicate or active-port termination"));
        }
        seen[port] = true;
        resolved[port] = load_impedance(
            &load.load,
            node_gas(net, exterior, media, port_nodes[port]).as_ref(),
            exterior,
            omega,
        )?;
    }
    if seen.iter().any(|&value| !value) {
        return Err(bad("network port termination omitted"));
    }

    let mut order: Vec<usize> = (0..net.nodes.len()).collect();
    order.sort_unstable_by(|&a, &b| net.nodes[a].id.cmp(&net.nodes[b].id));
    let mut ground = vec![false; net.nodes.len()];
    let mut node_loads = vec![None; net.nodes.len()];
    for &node in &order {
        if let NodeBoundary::Termination(load) = &net.nodes[node].boundary {
            node_loads[node] = load_impedance(
                load,
                node_gas(net, exterior, media, node).as_ref(),
                exterior,
                omega,
            )?;
            ground[node] = node_loads[node].is_some_and(|z| z == C64::ZERO);
        }
    }
    let mut slots = vec![None; net.nodes.len()];
    let mut terminal_slots = BTreeMap::new();
    let mut pressures = 0usize;
    for &node in &order {
        if let NodeBoundary::JunctionMultiport(feature) = &net.nodes[node].boundary {
            if feature.omega.to_bits() != omega.to_bits() {
                return Err(bad(
                    "junction multiport frequency differs from network query",
                ));
            }
            let mut terminals = feature.terminals.clone();
            terminals.sort();
            for terminal in terminals {
                terminal_slots.insert(terminal, pressures);
                pressures += 1;
            }
        } else if !ground[node] {
            slots[node] = Some(pressures);
            pressures += 1;
        }
    }
    let mut boundary_y = vec![C64::ZERO; pressures * pressures];
    for &node in &order {
        if let (Some(i), Some(z)) = (slots[node], node_loads[node]) {
            boundary_y[i * pressures + i] = z.recip();
        }
        if let NodeBoundary::JunctionMultiport(feature) = &net.nodes[node].boundary {
            let d = feature.terminals.len();
            for (i, ti) in feature.terminals.iter().enumerate() {
                for (j, tj) in feature.terminals.iter().enumerate() {
                    let a = terminal_slots[ti];
                    let b = terminal_slots[tj];
                    boundary_y[a * pressures + b] =
                        boundary_y[a * pressures + b] + feature.y[i * d + j];
                }
            }
        }
    }
    let mut edges: Vec<&NetworkEdge> = net.edges.iter().collect();
    edges.sort_unstable_by(|a, b| a.id.cmp(&b.id));
    let mut wave_edges = Vec::with_capacity(edges.len());
    for edge in edges {
        let interior = media.map_or(exterior, |binding| {
            binding
                .edge_gas(&edge.id)
                .expect("complete edge identities admitted")
        });
        let (two_port, condition, _) = network_chain_scattering(
            &edge.duct,
            interior,
            exterior,
            omega,
            loss,
            edge.corrections,
        )?;
        let endpoint = |node: usize, end: EdgeEnd| {
            if matches!(net.nodes[node].boundary, NodeBoundary::JunctionMultiport(_)) {
                terminal_slots
                    .get(&NetworkTerminal {
                        edge_id: edge.id.clone(),
                        end,
                    })
                    .copied()
            } else {
                slots[node]
            }
        };
        wave_edges.push(WaveEdge {
            ends: [
                endpoint(edge.from, EdgeEnd::From),
                endpoint(edge.to, EdgeEnd::To),
            ],
            two_port,
            condition,
        });
    }
    let port_slots: Vec<usize> = port_nodes
        .iter()
        .map(|&node| slots[node].expect("external port pressure is retained"))
        .collect();
    let (pressures, inflows) = loaded_wave_response(
        &boundary_y,
        pressures,
        &wave_edges,
        &labels,
        &port_slots,
        &conditions,
        &resolved,
        &coupled_rows,
        coupled,
    )?;
    Ok(LoadedFields {
        labels: labels.iter().map(|label| (*label).to_owned()).collect(),
        conditions,
        pressures,
        inflows,
    })
}

// Scalar/coupled loads are constitutive equations, not admittance stamps.
// Keeping their external currents permits a pressure release or a rank-deficient
// passive Z without an inverse, and closes the source/load problem at its pole.
#[allow(clippy::too_many_arguments)]
fn loaded_wave_response(
    boundary_y: &[C64],
    pressures: usize,
    edges: &[WaveEdge],
    labels: &[&str],
    port_slots: &[usize],
    conditions: &[Option<PortExcitation>],
    resolved: &[Option<C64>],
    coupled_rows: &[Option<usize>],
    coupled: Option<&CoupledTermination>,
) -> Result<(Vec<C64>, Vec<C64>), DuctError> {
    let source_flows: Vec<Option<C64>> = conditions
        .iter()
        .map(|condition| condition.and_then(PortExcitation::flow))
        .collect();
    // Pressure-source currents are unknowns just like external-load currents.
    // Keep their constitutive equations instead of dividing by Zs, including
    // an ideal pressure source with Zs=0 and rank-deficient passive loads.
    let mut retained_ports: Vec<usize> = (0..labels.len())
        .filter(|&port| source_flows[port].is_none())
        .collect();
    retained_ports.sort_unstable_by(|&a, &b| labels[a].cmp(labels[b]));
    let (mut a, n, reference) = wave_equations(boundary_y, pressures, edges, retained_ports.len())?;
    let mut current_slots = vec![None; labels.len()];
    let first_current = pressures + 2 * edges.len();
    for (offset, &port) in retained_ports.iter().enumerate() {
        current_slots[port] = Some(first_current + offset);
    }
    for &port in &retained_ports {
        let row = current_slots[port].expect("unprescribed port current is retained");
        a[port_slots[port] * n + row] = C64::from_re(-1.0);
        if let (Some(load), Some(load_row)) = (coupled, coupled_rows[port]) {
            a[row * n + port_slots[port]] = C64::ONE;
            for &other in &retained_ports {
                if let Some(load_col) = coupled_rows[other] {
                    let z = load.z[load_row * load.labels.len() + load_col];
                    let column = current_slots[other].expect("coupled current is retained");
                    a[row * n + column] = C64::new(z.re / reference, z.im / reference);
                }
            }
        } else if let Some(z) = resolved[port] {
            a[row * n + port_slots[port]] = C64::ONE;
            a[row * n + row] = C64::new(z.re / reference, z.im / reference);
        } else {
            // Explicitly closed load: U=0, with pressure unconstrained.
            a[row * n + row] = C64::ONE;
        }
    }
    let mut row_scales = vec![0.0; n];
    let mut pressure_scales = vec![0.0; n];
    for (row, values) in a.chunks_exact_mut(n).enumerate() {
        let scale = values.iter().map(|z| z.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale == 0.0 || values.iter().any(|&z| !finite(z)) {
            return Err(DuctError::Singular);
        }
        for value in values {
            *value = C64::new(value.re / scale, value.im / scale);
        }
        row_scales[row] = reference / scale;
        pressure_scales[row] = 1.0 / scale;
    }
    let (inv, _) = inverse(&a, n)?;
    let mut driven: Vec<usize> = (0..labels.len())
        .filter(|&port| conditions[port].is_some())
        .collect();
    driven.sort_unstable_by(|&a, &b| labels[a].cmp(labels[b]));
    let unit_source =
        (driven.len() == 1 && source_flows[driven[0]] == Some(C64::ONE)).then_some(driven[0]);
    let solution = |row: usize| {
        unit_source.map_or_else(
            || {
                coherent_sum(driven.iter().map(|&port| {
                    if let Some(flow) = source_flows[port] {
                        let slot = port_slots[port];
                        inv[row * n + slot] * flow.scale(row_scales[slot])
                    } else {
                        let slot =
                            current_slots[port].expect("pressure-source current is retained");
                        let (pressure, _) = conditions[port]
                            .and_then(PortExcitation::pressure)
                            .expect("driven pressure boundary");
                        inv[row * n + slot] * pressure.scale(pressure_scales[slot])
                    }
                }))
            },
            |active| {
                let slot = port_slots[active];
                inv[row * n + slot].scale(row_scales[slot])
            },
        )
    };
    let mut port_pressures = vec![C64::ZERO; labels.len()];
    let mut port_inflows: Vec<C64> = source_flows
        .iter()
        .map(|value| value.unwrap_or(C64::ZERO))
        .collect();
    for port in 0..labels.len() {
        port_pressures[port] = solution(port_slots[port]);
        if let Some(row) = current_slots[port] {
            let current = solution(row);
            port_inflows[port] = C64::new(current.re / reference, current.im / reference);
        }
    }
    if port_pressures
        .iter()
        .chain(&port_inflows)
        .any(|&z| !finite(z))
    {
        return Err(DuctError::Singular);
    }
    Ok((port_pressures, port_inflows))
}
