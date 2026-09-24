//! Geometry-bound graph input. Every node lowers to an existing physical owner;
//! no branch is replaced with a second independent voice or a fitted reflection.
use super::*;
use crate::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use crate::bernoulli_aperture::wall::{WallPatch, WallPin};
use fs_vfit::impedance::SeriesImpedanceSpec;
use crate::bernoulli_aperture::viscothermal::{
    ViscothermalSectionSpec, ViscothermalSection, with_viscothermal_sections,
};
use fs_vfit::relaxation::RelaxationImpedanceSpec;

const MAX_DUCT_NODES: usize = 64;
const MAX_DUCT_SECTIONS: usize = 128;
const MAX_WAVE_BYTES: usize = 64 * 1024 * 1024;

pub(super) enum DuctInput {
    // Keep the original source's arithmetic and runtime when no graph is named.
    Tube { spec: UniformTubeSpec, radiation_band: Option<f64>, loss: Option<ViscothermalSectionSpec> },
    Graph { nodes: Vec<NodeInput>, sections: Vec<TubeSection>, max_bytes: usize, losses: Vec<Option<ViscothermalSectionSpec>> },
}

pub(super) enum NodeInput {
    Existing(NetworkNode),
    Radiation(f64),
    Cavity(HelmholtzLoadSpec),
    Wall(WallPatch),
}

pub(super) struct PreparedDuct {
    tube: Option<UniformTubeSpec>,
    graph: Option<TubeNetworkSpec>,
    pub observation: ApertureObservation,
    pub radiation: Vec<(usize, BaffledRadiationLoad)>,
    pub losses: Vec<ViscothermalSection>,
}

impl DuctInput {
    pub fn read(r: &mut Reader<'_>, air: &GasState, gate: &CancelGate)
        -> Result<Self, PlateValveInputError>
    {
        let mut row = r.next()?;
        match row.word()? {
            "tube" => {
                let length_m = row.scalar()?;
                let radius_m = row.scalar()?;
                let terminal = row.word()?;
                let (terminal_reflection, radiation_band) = if terminal == "baffled-low-ka" {
                    // Geometry helper only; prepare replaces this before binding.
                    (0.0, Some(row.scalar()?))
                } else {
                    let reflection: f64 = terminal.parse()
                        .map_err(|_| bad(r.line, "tube requires a numeric reflectance or baffled-low-ka BAND_HZ"))?;
                    if !reflection.is_finite() || reflection.abs() > 1.0 {
                        return Err(bad(r.line, "tube reflectance must be finite and passive"));
                    }
                    (reflection, None)
                };
                let spec = UniformTubeSpec { length_m, radius_m, terminal_reflection,
                    max_length_error_m: row.scalar()?, max_wave_memory_bytes: row.count(MAX_WAVE_BYTES)?,
                    sound_speed_m_s: air.sound_speed };
                let loss = read_loss(&mut row)?;
                row.finish()?;
                Ok(Self::Tube { spec, radiation_band, loss })
            }
            "network" => {
                let count = row.count(MAX_DUCT_NODES)?;
                let section_count = row.count(MAX_DUCT_SECTIONS)?;
                let max_bytes = row.count(MAX_WAVE_BYTES)?;
                row.finish()?;
                if count < 2 || section_count == 0 {
                    return Err(bad(r.line, "network requires 2..=64 nodes and 1..=128 physical sections"));
                }
                let mut nodes = Vec::with_capacity(count);
                for _ in 0..count {
                    checkpoint(gate)?;
                    let mut row = r.row("duct_node")?;
                    let input = match row.word()? {
                        "inlet" => NodeInput::Existing(NetworkNode::Inlet),
                        "junction" => NodeInput::Existing(NetworkNode::Junction),
                        "reflection" => {
                            let reflection = row.scalar()?;
                            if reflection.abs() > 1.0 { return Err(bad(r.line, "terminal reflectance must be passive")); }
                            NodeInput::Existing(NetworkNode::Termination { reflection })
                        }
                        "baffled-low-ka" => NodeInput::Radiation(row.scalar()?),
                        "cavity" => NodeInput::Cavity(HelmholtzLoadSpec {
                            volume_m3: row.scalar()?, neck_radius_m: row.scalar()?,
                            effective_neck_length_m: row.scalar()?, resistance_pa_s_m3: row.scalar()?,
                        }),
                        "wall" => NodeInput::Wall(WallPatch {
                            area_m2: row.scalar()?, wall: WallPin {
                                surface_density: row.scalar()?, stiffness_per_area: row.scalar()?, resistance: row.scalar()?,
                            },
                        }),
                        kind @ ("impedance" | "series" | "shunt") => {
                            let load = series_load(&mut row)?;
                            NodeInput::Existing(match kind {
                                "impedance" => NetworkNode::Impedance { load },
                                "series" => NetworkNode::Series { load: relaxation(load)? },
                                _ => NetworkNode::Shunt { load: relaxation(load)? },
                            })
                        }
                        _ => return Err(bad(r.line, "unsupported physical duct node")),
                    };
                    row.finish()?;
                    nodes.push(input);
                }
                let mut sections = Vec::with_capacity(section_count);
                let mut losses = Vec::with_capacity(section_count);
                for _ in 0..section_count {
                    checkpoint(gate)?;
                    let mut row = r.row("duct_section")?;
                    let nodes = [row.parse()?, row.parse()?];
                    if nodes[0] == nodes[1] || nodes.iter().any(|&i| i >= count) {
                        return Err(bad(r.line, "duct section needs distinct existing endpoints"));
                    }
                    let section = TubeSection { nodes, length_m: row.scalar()?, radius_m: row.scalar()?,
                        max_length_error_m: row.scalar()? };
                    losses.push(read_loss(&mut row)?);
                    row.finish()?;
                    if section.length_m <= 0.0 || section.radius_m <= 0.0 || section.max_length_error_m < 0.0 {
                        return Err(bad(r.line, "duct sections require positive dimensions and nonnegative length error"));
                    }
                    sections.push(section);
                }
                Ok(Self::Graph { nodes, sections, max_bytes, losses })
            }
            _ => Err(bad(r.line, "expected an explicit tube or network after ambient")),
        }
    }

    pub fn prepare(self, air: &GasState, dt: f64, observation: ApertureObservation, gate: &CancelGate)
        -> Result<PreparedDuct, PlateValveInputError>
    {
        let physics = PlateValveInputError::Physics;
        match self {
            Self::Tube { spec, radiation_band, loss } => {
                if !matches!(observation, ApertureObservation::Inlet | ApertureObservation::TubeTerminal
                    | ApertureObservation::TubeBaffled(_)) {
                    return Err(bad(0, "a tube observation must name inlet, terminal or baffled-outlet"));
                }
                if radiation_band.is_none() && loss.is_none() {
                    return Ok(PreparedDuct { tube: Some(spec), graph: None, observation, radiation: vec![], losses: vec![] });
                }
                let mut radiation = Vec::new();
                let terminal = if let Some(band) = radiation_band {
                    let load = BaffledRadiationLoad::new(spec.radius_m, air.density, air.sound_speed, dt, band, gate)
                        .map_err(physics)?;
                    radiation.push((1, load));
                    load.termination()
                } else { NetworkNode::Termination { reflection: spec.terminal_reflection } };
                let observation = match observation {
                    ApertureObservation::Inlet => observation,
                    ApertureObservation::TubeTerminal => ApertureObservation::NetworkNode(1),
                    ApertureObservation::TubeBaffled(receiver) => {
                        if let Some(band) = radiation_band {check_receiver_band(receiver, band)?;}
                        ApertureObservation::NetworkBaffled { node: 1, receiver }
                    }
                    _ => unreachable!("tube observation was admitted"),
                };
                let mut prepared = PreparedDuct { tube: None, graph: Some(TubeNetworkSpec {
                    nodes: vec![NetworkNode::Inlet, terminal],
                    sections: vec![TubeSection { nodes: [0, 1], length_m: spec.length_m,
                        radius_m: spec.radius_m, max_length_error_m: spec.max_length_error_m }],
                    sound_speed_m_s: air.sound_speed, max_wave_memory_bytes: spec.max_wave_memory_bytes,
                }), observation, radiation, losses: vec![] };
                prepared.install_losses(air, dt, &[loss], gate)?;
                Ok(prepared)
            }
            Self::Graph { nodes, sections, max_bytes, losses } => {
                match observation {
                    ApertureObservation::Inlet => {},
                    ApertureObservation::NetworkNode(node) if node < nodes.len() => {},
                    ApertureObservation::NetworkBaffled { node, .. } if node < nodes.len() => {
                        if !matches!(&nodes[node], NodeInput::Radiation(_) | NodeInput::Existing(
                            NetworkNode::Termination { .. } | NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. })) {
                            return Err(bad(0, "exterior graph observation must name an outlet, not an enclosed cavity or interior node"));
                        }
                    }
                    _ => return Err(bad(0, "a network needs inlet, network-node INDEX or network-baffled INDEX and receiver")),
                }
                let mut lowered = Vec::with_capacity(nodes.len());
                let mut radiation = Vec::new();
                for (node, input) in nodes.into_iter().enumerate() {
                    checkpoint(gate)?;
                    lowered.push(match input {
                        NodeInput::Existing(kind) => kind,
                        NodeInput::Cavity(cavity) => cavity.termination(air.density, air.sound_speed).map_err(physics)?,
                        NodeInput::Wall(wall) => wall.shunt().map_err(physics)?,
                        NodeInput::Radiation(band) => {
                            let mut incident = sections.iter().filter(|s| s.nodes.contains(&node));
                            let section = incident.next().ok_or_else(|| bad(0, "radiating terminal is disconnected"))?;
                            if incident.next().is_some() { return Err(bad(0, "radiating terminal must meet exactly one physical section")); }
                            let load = BaffledRadiationLoad::new(section.radius_m, air.density, air.sound_speed, dt, band, gate)
                                .map_err(physics)?;
                            if let ApertureObservation::NetworkBaffled { node: selected, receiver } = observation {
                                if node == selected { check_receiver_band(receiver, band)?; }
                            }
                            radiation.push((node, load));
                            load.termination()
                        }
                    });
                }
                let graph = TubeNetworkSpec { nodes: lowered, sections, sound_speed_m_s: air.sound_speed,
                    max_wave_memory_bytes: max_bytes };
                // Admit the actual pressure area before deriving the moving plate.
                // Full graph topology and storage stay owned by ApertureNetwork.
                graph.inlet_impedance(air.density).map_err(physics)?;
                let mut prepared = PreparedDuct { tube: None, graph: Some(graph), observation, radiation, losses: vec![] };
                prepared.install_losses(air, dt, &losses, gate)?;
                Ok(prepared)
            }
        }
    }
}

impl PreparedDuct {
    fn install_losses(&mut self, air: &GasState, dt: f64,
        assignments: &[Option<ViscothermalSectionSpec>], gate: &CancelGate) -> Result<(), PlateValveInputError>
    {
        if assignments.iter().all(Option::is_none) {return Ok(());}
        if let ApertureObservation::NetworkBaffled {receiver, ..} = self.observation {
            for options in assignments.iter().flatten() {
                if receiver.maximum_frequency_hz > options.maximum_frequency_hz {
                    return Err(bad(0, "receiver band cannot exceed any selected viscothermal section band"));
                }
            }
        }
        let graph = self.graph.as_ref().expect("loss selection uses the existing network owner");
        let selections: Vec<_> = assignments.iter().enumerate().filter_map(|(section, selected)|
            selected.map(|s| ViscothermalSectionSpec { section, ..s })).collect();
        let (lowered, reports) = with_viscothermal_sections(graph.clone(), air, dt, &selections, gate)
            .map_err(PlateValveInputError::Physics)?;
        self.graph = Some(lowered);
        self.losses = reports;
        Ok(())
    }

    pub fn inlet_impedance(&self, density: f64) -> Result<f64, PlateValveInputError> {
        match (&self.tube, &self.graph) {
            (Some(tube), None) => tube.characteristic_impedance(density),
            (None, Some(graph)) => graph.inlet_impedance(density),
            _ => unreachable!("exactly one propagation owner is selected"),
        }.map_err(PlateValveInputError::Physics)
    }

    // The convenience field describes the observed outlet. Internal observations
    // retain it only for a unique radiating terminal; the full list is separate.
    pub fn observed_radiation(&self) -> Option<BaffledRadiationLoad> {
        if let ApertureObservation::NetworkBaffled { node, .. } = self.observation {
            self.radiation.iter().find(|(i, _)| *i == node).map(|(_, load)| *load)
        } else if self.radiation.len() == 1 { Some(self.radiation[0].1) } else { None }
    }

    pub fn bind(self, valve: DynamicAperture) -> Result<CoupledAperture, PlateValveInputError> {
        match (self.tube, self.graph) {
            (Some(tube), None) => ApertureTube::new(valve, tube).map(CoupledAperture::Tube),
            (None, Some(graph)) => ApertureNetwork::new(valve, graph).map(CoupledAperture::Network),
            _ => unreachable!("exactly one propagation owner is selected"),
        }.map_err(PlateValveInputError::Physics)
    }
}

fn series_load(row: &mut Row<'_>) -> Result<SeriesImpedanceSpec, PlateValveInputError> {
    let resistance_pa_s_m3 = row.scalar()?;
    let inertance_pa_s2_m3 = row.scalar()?;
    let value = row.word()?;
    let compliance_m3_pa = if value == "none" { None } else {
        Some(value.parse::<f64>().map_err(|_| bad(row.line, "compliance needs positive cubic metres/Pa or none"))?)
    };
    let load = SeriesImpedanceSpec { resistance_pa_s_m3, inertance_pa_s2_m3, compliance_m3_pa };
    load.validate().map_err(|_| bad(row.line, "duct impedance requires finite passive R/L and positive optional C"))?;
    Ok(load)
}
fn relaxation(base: SeriesImpedanceSpec) -> Result<RelaxationImpedanceSpec, PlateValveInputError> {
    RelaxationImpedanceSpec::new(base, &[])
        .map_err(|_| bad(0, "interior duct impedance was refused"))
}
fn check_receiver_band(receiver: CircularOutletReceiver, band: f64) -> Result<(), PlateValveInputError> {
    if receiver.maximum_frequency_hz > band {
        return Err(bad(0, "exterior receiver band cannot exceed the selected radiation-load band"));
    }
    Ok(())
}

// Optional suffix on the original tube/section record. It names numerical
// approximation choices only: all transport coefficients come from `ambient`.
fn read_loss(row: &mut Row<'_>) -> Result<Option<ViscothermalSectionSpec>, PlateValveInputError> {
    let Some(kind) = row.fields.next() else {return Ok(None);};
    if kind != "viscothermal" {return Err(bad(row.line, "expected viscothermal MIN_HZ MAX_HZ CELLS ARMS or end of section"));}
    let selected = ViscothermalSectionSpec {
        section: 0, // Assigned from the original section address at preparation.
        minimum_frequency_hz: row.scalar()?, maximum_frequency_hz: row.scalar()?,
        cells: row.count(128)?,
    };
    // The shared owner has an explicit eight-arm spectrum. Do not silently
    // accept and ignore a caller's different finite-order request.
    if row.count(8)? != 8 {return Err(bad(row.line, "the shared viscothermal owner requires exactly eight arms"));}
    Ok(Some(selected))
}
