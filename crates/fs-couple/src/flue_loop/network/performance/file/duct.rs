//! Geometry and gas data lowered into the existing shared tube-network owner.
//! No propagation, scattering, load or material formula is repeated here.
use super::*;
use crate::flue_loop::network::performance::FlueRadiatingOutlet;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::RelaxationImpedanceSpec;

pub(super) struct Input {
    pub graph: TubeNetworkSpec,
    pub gases: Vec<GasState>,
    pub inner: Vec<HoleInnerLength>,
    explicit_network: bool,
}

/// An open side hole's inner length, from `inner BORE_RADIUS_M` on its
/// chimney's `duct_section` row (before any `gas`):
/// `t_i = b (0.82 - 0.193 d - 1.09 d^2 + 1.27 d^3 - 0.71 d^4)`, `d = b / a`,
/// with b the chimney's radius and a the main bore's radius at the hole.
/// This is Dubos et al. (1999) as corrected by Dalmont et al. (2002), in the
/// form Garcia Mayen et al. give it (HAL hal-03230052, 2021, eq. 11) and
/// Ernoult and Fabre (2017) use (eq. 17: twice it is their right-angle bend).
/// It is the mass of the junction's evanescent modes, in the hole's area.
///
/// Reference plane: the chimney row's length is the hole's height on its axis
/// from the main bore's wall (a card's wall thickness at the hole), and the
/// bore's sections run through the junction at the hole's axis. That is the
/// plane of Garcia Mayen et al.'s Fig. 1, whose shunt length is
/// `t_s = t_i + t + t_m + t_r` (eq. 14). This row adds `t_i` once, to the
/// chimney's length. It does not add the matching volume `t_m` (their eq. 12),
/// the series length `t_a` (eq. 10), or anything at the outer end, whose load
/// is the outlet's own. The row needs a side hole: a chimney from a junction
/// of at least three sections to an open end (reflection -1), narrower than
/// the bore (`0 < d <= 1`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HoleInnerLength {
    /// The chimney's section index.
    pub section: usize,
    /// b, the hole's radius [m].
    pub hole_radius_m: f64,
    /// a, the main bore's radius at the hole, from the row [m].
    pub bore_radius_m: f64,
    /// The chimney as the row gives it: the height on the hole's axis [m].
    pub chimney_m: f64,
    /// `t_i` [m], added to the chimney's length in the graph.
    pub inner_length_m: f64,
}

impl HoleInnerLength {
    /// `t_i` [m] for a hole of radius b on a bore of radius a; `None` unless
    /// `0 < b / a <= 1` with both finite.
    #[must_use]
    pub fn law(hole_radius_m: f64, bore_radius_m: f64) -> Option<f64> {
        let d = hole_radius_m / bore_radius_m;
        (hole_radius_m.is_finite() && hole_radius_m > 0.0 && bore_radius_m.is_finite()
            && d > 0.0 && d <= 1.0)
            .then(|| hole_radius_m * ((((-0.71 * d + 1.27) * d - 1.09) * d - 0.193) * d + 0.82))
    }
}

impl Input {
    pub fn read(rows: &mut Rows<'_>, default: GasState,
        sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>,
        gate: &CancelGate) -> Result<Self, FlueInputError>
    {
        let explicit_network = rows.next_is("network");
        let mut inner = Vec::new();
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
            let mut requests = Vec::new();
            for index in 0..section_count {
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
                if row.fields.clone().next() == Some("inner") {
                    row.fields.next();
                    requests.push((index, row.scalar()?, row.line));
                }
                let gas = if row.fields.clone().next() == Some("gas") {
                    row.fields.next();
                    row.gas(sources, records)?
                } else { default };
                row.end()?;
                sections.push(section);
                gases.push(gas);
            }
            for (index, bore_radius_m, line) in requests {
                let TubeSection { nodes: [near, far], length_m, radius_m, .. } = sections[index];
                let open = matches!(nodes[far],
                    NetworkNode::Termination { reflection } if reflection.to_bits() == (-1f64).to_bits());
                let degree = sections.iter().filter(|s| s.nodes.contains(&near)).count();
                if !open || !matches!(nodes[near], NetworkNode::Junction) || degree < 3 {
                    return Err(bad(line, "duct_section inner needs a side hole: a chimney from a junction of at least three sections to an open end (reflection -1)"));
                }
                let Some(inner_length_m) = HoleInnerLength::law(radius_m, bore_radius_m) else {
                    return Err(bad(line, "duct_section inner needs a hole no wider than the bore: 0 < hole radius / bore radius <= 1"));
                };
                sections[index].length_m = length_m + inner_length_m;
                inner.push(HoleInnerLength { section: index, hole_radius_m: radius_m, bore_radius_m,
                    chimney_m: length_m, inner_length_m });
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
        Ok(Self { graph, gases, inner, explicit_network })
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
