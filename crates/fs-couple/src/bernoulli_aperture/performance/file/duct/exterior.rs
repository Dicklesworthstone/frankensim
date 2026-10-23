//! One resolved exterior for the existing compact loads and flow observation.
//! Original gas, loss, cavity, plate and control owners remain unchanged.
use super::*;

impl DuctInput {
    /// Lower the original duct first, reserving declared radiating endpoints for
    /// the selected exterior. Never build a false local-gas radiation trial: its
    /// ka refusal need not apply to the actual exterior. All other loads remain
    /// authored. Geometry and regional loss subdivision use the original owner.
    pub fn prepare_in_exterior(self, air: &GasState, dt: f64,
        observation: ApertureObservation, exterior: Option<&GasState>, gate: &CancelGate)
        -> Result<PreparedDuct, PlateValveInputError>
    {
        let Some(exterior) = exterior else { return self.prepare(air, dt, observation, gate); };
        checkpoint(gate)?;
        exterior.validate_coherence().map_err(|source| PlateValveInputError::Gas { line: 0, source })?;
        if exterior.pressure != air.pressure {
            return Err(bad(0, "exterior and ambient must share static pressure; mean interface flow is not modeled"));
        }
        // The legacy radiation tube already lowers to this same network. Only
        // the typed regional binding is made explicit; the gas is INTERNAL air.
        let (mut source, observation) = match self {
            Self::Dg(_) => return Err(bad(0, "dg_tube has no bound exterior radiation operator; its inlet pressure cannot be an exterior receiver")),
            Self::Tube { spec, radiation_band: Some(band), loss } => {
                let observation = match observation {
                    ApertureObservation::Inlet => observation,
                    ApertureObservation::TubeTerminal => ApertureObservation::NetworkNode(1),
                    ApertureObservation::TubeBaffled(receiver) => ApertureObservation::NetworkBaffled { node: 1, receiver },
                    _ => return Err(bad(0, "a tube needs its original inlet, terminal or baffled-outlet observation")),
                };
                (Self::Graph {
                    nodes: vec![NodeInput::Existing(NetworkNode::Inlet), NodeInput::Radiation(band)],
                    sections: vec![TubeSection { nodes: [0, 1], length_m: spec.length_m,
                        radius_m: spec.radius_m, max_length_error_m: spec.max_length_error_m }],
                    max_bytes: spec.max_wave_memory_bytes, losses: vec![loss], gases: Some(vec![*air]),
                }, observation)
            }
            Self::Tube { .. } => return Err(bad(0, "explicit exterior requires a declared baffled-low-ka radiation load; numeric reflection is not replaced")),
            graph @ Self::Graph { .. } => (graph, observation),
        };
        let Self::Graph { nodes, sections, gases, .. } = &mut source else { unreachable!() };
        let mut radiating = Vec::new();
        for (node, input) in nodes.iter_mut().enumerate() {
            if let NodeInput::Radiation(band) = input {
                radiating.push((node, *band));
                *input = NodeInput::Existing(NetworkNode::Termination { reflection: 0.0 });
            }
        }
        if radiating.is_empty() {
            return Err(bad(0, "explicit exterior requires at least one declared baffled-low-ka radiation endpoint"));
        }
        if let ApertureObservation::NetworkBaffled { node, receiver } = observation {
            let band = radiating.iter().find(|(id, _)| *id == node)
                .ok_or_else(|| bad(0, "explicit exterior observation must select a declared radiation endpoint, not an arbitrary impedance"))?.1;
            check_receiver_band(receiver, band)?;
        }
        let internal = gases.get_or_insert_with(|| vec![*air; sections.len()]);
        if internal.iter().any(|gas| gas.pressure != exterior.pressure) {
            return Err(bad(0, "every internal region must share the exterior static pressure"));
        }
        let mut prepared = source.prepare(air, dt, observation, gate)?;
        let graph = prepared.graph.as_mut().expect("exterior source uses the network owner");
        for (node, band) in radiating {
            checkpoint(gate)?;
            let (section, _) = terminal_gas(&graph.sections, prepared.gases.as_deref(), node, air)?;
            let load = BaffledRadiationLoad::new(section.radius_m, exterior.density,
                exterior.sound_speed, dt, band, gate).map_err(PlateValveInputError::Physics)?;
            graph.nodes[node] = load.termination();
            prepared.radiation.push((node, load));
        }
        checkpoint(gate)?;
        Ok(prepared)
    }
}

#[cfg(test)]
mod tests;
