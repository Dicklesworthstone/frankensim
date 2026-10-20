//! Geometry and gas data lowered into the existing shared tube-network owner.
//! No propagation, scattering, load or material formula is repeated here.
use super::*;
use crate::flue_loop::network::performance::FlueRadiatingOutlet;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::RelaxationImpedanceSpec;

pub(super) struct Input {
    pub graph: TubeNetworkSpec,
    pub gases: Vec<GasState>,
    explicit_network: bool,
}

impl Input {
    pub fn read(rows: &mut Rows<'_>, default: GasState,
        sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>,
        gate: &CancelGate) -> Result<Self, FlueInputError>
    {
        let explicit_network = rows.next_is("network");
        let (nodes, sections, gases) = if explicit_network {
            let mut row = rows.row("network")?;
            let node_count: usize = row.number()?;
            let section_count: usize = row.number()?;
            row.end()?;
            if !(2..=64).contains(&node_count) || !(1..=128).contains(&section_count) {
                return Err(bad(rows.line, "network requires 2..=64 nodes and 1..=128 sections"));
            }
            let mut nodes = Vec::with_capacity(node_count);
            for _ in 0..node_count {
                poll(gate)?;
                let mut row = rows.row("duct_node")?;
                let node = match row.fields.next() {
                    Some("inlet") => NetworkNode::Inlet,
                    Some("driven") => NetworkNode::Driven,
                    Some("junction") => NetworkNode::Junction,
                    Some("reflection") => {
                        let reflection = row.scalar()?;
                        if reflection.abs() > 1.0 {
                            return Err(bad(row.line, "terminal reflection must be passive"));
                        }
                        NetworkNode::Termination { reflection }
                    }
                    Some(kind @ ("impedance" | "series" | "shunt")) => {
                        let resistance = row.scalar()?;
                        let inertance = row.scalar()?;
                        let compliance = row.scalar()?;
                        // Zero explicitly denotes absent compliance; it is not
                        // an inferred open/closed boundary or a numerical fit.
                        let load = SeriesImpedanceSpec { resistance_pa_s_m3: resistance,
                            inertance_pa_s2_m3: inertance,
                            compliance_m3_pa: (compliance != 0.0).then_some(compliance) };
                        load.validate().map_err(|_| bad(row.line, "invalid passive R-L-C duct load"))?;
                        match kind {
                            "impedance" => NetworkNode::Impedance { load },
                            kind => {
                                let load = RelaxationImpedanceSpec::new(load, &[])
                                    .map_err(|_| bad(row.line, "invalid passive interior duct load"))?;
                                if kind == "series" { NetworkNode::Series { load } }
                                else { NetworkNode::Shunt { load } }
                            }
                        }
                    }
                    _ => return Err(bad(row.line, "unsupported physical duct node")),
                };
                row.end()?;
                nodes.push(node);
            }
            let mut sections = Vec::with_capacity(section_count);
            let mut gases = Vec::with_capacity(section_count);
            for _ in 0..section_count {
                poll(gate)?;
                let mut row = rows.row("duct_section")?;
                let section = TubeSection { nodes: [row.number()?, row.number()?],
                    length_m: row.scalar()?, radius_m: row.scalar()?, max_length_error_m: row.scalar()? };
                if section.nodes[0] == section.nodes[1]
                    || section.nodes.iter().any(|&node| node >= node_count)
                    || section.length_m <= 0.0 || section.radius_m <= 0.0
                    || section.max_length_error_m < 0.0 {
                    return Err(bad(row.line, "duct sections require distinct existing endpoints and positive geometry"));
                }
                let gas = if row.fields.clone().next() == Some("gas") {
                    row.fields.next();
                    row.gas(sources, records)?
                } else { default };
                row.end()?;
                sections.push(section);
                gases.push(gas);
            }
            (nodes, sections, gases)
        } else {
            let mut row = rows.row("tube")?;
            let length_m = row.scalar()?;
            let radius_m = row.scalar()?;
            let max_length_error_m = row.scalar()?;
            let reflection = row.scalar()?;
            row.end()?;
            (vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
                vec![TubeSection { nodes: [0, 1], length_m, radius_m, max_length_error_m }],
                vec![default])
        };
        let mut graph = TubeNetworkSpec { nodes, sections,
            sound_speed_m_s: default.sound_speed, max_wave_memory_bytes: 0 };
        let inlet = graph.inlet_section_index().map_err(|source|
            FlueInputError::Physics(FlueError::Physics(source)))?;
        graph.sound_speed_m_s = gases[inlet].sound_speed;
        graph.validate_section_gases(&gases).map_err(|source|
            FlueInputError::Physics(FlueError::Physics(source)))?;
        Ok(Self { graph, gases, explicit_network })
    }

    pub fn read_losses(&self, rows: &mut Rows<'_>) -> Result<Vec<ViscothermalSectionSpec>, FlueInputError> {
        let record = rows.record("viscothermal")?;
        if record == "none" { return Ok(Vec::new()); }
        let mut row = Row { fields: record.split_ascii_whitespace(), line: rows.line };
        if row.fields.clone().next() == Some("sections") {
            row.fields.next();
            let count: usize = row.number()?;
            row.end()?;
            if count > self.graph.sections.len() {
                return Err(bad(rows.line, "viscothermal selection exceeds the physical section count"));
            }
            let mut selections = Vec::<ViscothermalSectionSpec>::with_capacity(count);
            for _ in 0..count {
                let mut row = rows.row("loss_section")?;
                let section = row.number()?;
                if section >= self.graph.sections.len() || selections.iter().any(|s| s.section == section) {
                    return Err(bad(row.line, "loss_section must name one distinct existing section"));
                }
                let loss = ViscothermalSectionSpec { section, cells: row.number()?,
                    minimum_frequency_hz: row.scalar()?, maximum_frequency_hz: row.scalar()? };
                row.end()?;
                selections.push(loss);
            }
            Ok(selections)
        } else {
            if self.explicit_network {
                return Err(bad(row.line, "network losses require none or sections COUNT with addressed loss_section records"));
            }
            let loss = ViscothermalSectionSpec { section: 0, cells: row.number()?,
                minimum_frequency_hz: row.scalar()?, maximum_frequency_hz: row.scalar()? };
            row.end()?;
            Ok(vec![loss])
        }
    }

    pub fn read_outlets(&self, rows: &mut Rows<'_>) -> Result<(usize, Vec<FlueRadiatingOutlet>), FlueInputError> {
        if !rows.next_is("outlets") && !self.explicit_network { return Ok((1, Vec::new())); }
        let mut row = rows.row("outlets")?;
        let primary: usize = row.number()?;
        let additional: usize = row.number()?;
        row.end()?;
        if primary >= self.graph.nodes.len()
            || additional >= crate::pcm_wav::baffled::MAX_RAYLEIGH_INPUTS {
            return Err(bad(rows.line, "outlets exceed the physical node or receiver input budget"));
        }
        let mut outlets = Vec::with_capacity(additional);
        for _ in 0..additional {
            let mut row = rows.row("outlet")?;
            let node: usize = row.number()?;
            if node >= self.graph.nodes.len() || node == primary
                || outlets.iter().any(|p: &FlueRadiatingOutlet| p.node == node) {
                return Err(bad(row.line, "outlet must name one distinct existing node"));
            }
            let center_m = [row.scalar()?, row.scalar()?];
            row.end()?;
            outlets.push(FlueRadiatingOutlet { node, center_m });
        }
        // The shared DuctPressure owner validates terminal degree, overlap,
        // common-pressure exterior and the geometry/band of every aperture.
        Ok((primary, outlets))
    }
}
