//! The nasal tract and the maxillary sinuses beside the card's oral tract (N3; owner order MB 15368; design 17372 as
//! corrected by 17432 and 17442; rulings MB 17768). All rows are optional, and a card without them never reaches this
//! file.
//!
//! The rows carry measured geometry inline, as `section` rows carry the oral chart; no data file is read here:
//!
//! ```text
//! nasal_source <identifier without spaces>          (needed with any nasal row: the nasal tables)
//! nasal_licence <free text to end of line>
//! velum <after_section> <port_area_m2> <attach_source without spaces>
//! nasal_section <area_m2> <length_m>                 (repeat, velum first: the common passage)
//! nasal_pair <left_area_m2> <right_area_m2> <length_m>   (repeat, choana first: the two passages)
//! sinus <left|right> <from_nostril_m> <ostium_area_m2> <ostium_length_m> <volume_m3>   (optional, one per side)
//! nostrils <x_m> <y_m> <source without spaces>       (the nostrils' place beside the mouth on the baffle)
//! ```
//!
//! The graph, every edge a chain of the sections' own cylinders (radius `sqrt(A / pi)`, as the oral chart's duct):
//! the oral tract split at the junction after section `after_section` (0-based); from that junction the nasal
//! passage, whose first section is the velopharyngeal port at `port_area_m2` (the nasal table's posterior section at
//! the port's area, MB 17768 ruling 2: the row carries no port length); at its end the choana junction, from which the
//! left and right passages run to the nostrils. A sinus opens off its passage `from_nostril_m` before the nostril, as
//! its ostium (a cylinder of the measured area and length) ending in the measured cavity's lumped compliance
//! ([`fs_duct::network::NetworkLoad::Cavity`]: the source's own model, no cavity shape; it holds below the cavity's
//! first internal mode). A passage section the ostium falls inside is split into two cylinders of the same area, which
//! is exact. Each nostril radiates at its passage's last radius with the card's own `termination`, as the lips do.
//! Every edge carries the card's `tract_wall` pin, and the tract's gas fills every edge.
//!
//! The island's line is the network's glottal line on the island's own ladder, against the island's own `rho c / S`
//! ([`GlottalIsland::new_with_tract_line`]); the lips' and each nostril's outward flow per incident wave are the same
//! solve's port flows ([`NetworkDrivingPoint::realize_port_flow_ir_regional`]) on the series mouth line's default
//! ladder's rungs and tolerances up to [`NASAL_FLOW_MAX_TAPS`], so a flow the series cap holds stops on the series
//! rung, and an open port's flows, whose slow tails still miss the tolerances at 65536 taps, get two more doublings.
//! The receiver takes the three openings on one baffle: the mouth at the origin and the nostrils at `(x, y)` and
//! `(-x, y)`, the baffle's `x` the subject's left. It is one plane: a nostril in front of the mouth's plane (the
//! nostrils row's source states by how much) is not modeled in front of it, and the openings do not load one another.
//!
//! Refused with a nasal row, by name: `tract_move` (the moving tract's lines are the series duct's), `head_sphere`,
//! `head_load` and `torso_sphere` (that load and that sphere are one mouth's), `tract_gas` (its table is one mouth's),
//! `fricative_noise` (its junction lines are the series duct's) and `receiver free_field_monopole` (one mouth's
//! monopole).

use super::*;
use crate::driving_point::network::NetworkDrivingPoint;
use crate::driving_point::outlet::MouthFlowRealization;
use crate::glottis::glottal_port_impedance;
use fs_duct::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkMedia, NetworkNode, NodeBoundary,
    PortTermination,
};
use fs_duct::{Duct, Segment};
use fs_ports::PortKind;

/// Every nasal row's key.
pub(super) const NASAL_KEYS: [&str; 7] =
    ["nasal_source", "nasal_licence", "velum", "nasal_section", "nasal_pair", "sinus", "nostrils"];

/// An ostium cut this close to a passage section's end is at that end [m]: rounding in `from_nostril_m` leaves no
/// sliver cylinder, and a sinus this close to either end of its passage is refused.
const CUT_SNAP_M: f64 = 1.0e-12;

/// Cap of the lips' and nostrils' flow ladders [taps]: the series line's rungs and tolerances, two doublings past its
/// cap. An open port's flows carry a slow, smooth tail that falls in late time (no content before time 0, no line):
/// on the walled and the rigid test fixtures at 96 kHz the lips stop at 131072 taps and the left nostril still misses
/// the tap-change tolerance there. Each flow is a direct FIR, so every doubling doubles its cost per sample; past this
/// cap the ladder refuses by name.
const NASAL_FLOW_MAX_TAPS: usize = 262_144;

/// One `nasal_pair` row: the left and right passages' areas over one length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NasalPair {
    /// Left passage area [m^2].
    pub left_area_m2: f64,
    /// Right passage area [m^2].
    pub right_area_m2: f64,
    /// Length of both [m].
    pub length_m: f64,
}

/// One `sinus` row: a maxillary sinus off its side's passage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NasalSinus {
    /// Distance from the nostril back along the passage to the ostium [m].
    pub from_nostril_m: f64,
    /// Ostium area [m^2].
    pub ostium_area_m2: f64,
    /// Ostium length [m].
    pub ostium_length_m: f64,
    /// Cavity volume [m^3].
    pub volume_m3: f64,
}

/// The admitted nasal rows (module doc).
#[derive(Debug, Clone, PartialEq)]
pub struct NasalCard {
    /// `nasal_source`.
    pub source: String,
    /// `nasal_licence`.
    pub licence: String,
    /// The oral section (0-based) after which the port opens.
    pub after_section: usize,
    /// Velopharyngeal port area [m^2]: the area of the nasal passage's first section.
    pub port_area_m2: f64,
    /// Where on the card's oral tract the port opens: the source of `after_section`.
    pub attach_source: String,
    /// The common passage, velum first; section 0's own area is replaced by `port_area_m2`.
    pub sections: Vec<TractSection>,
    /// The two passages, choana first.
    pub pairs: Vec<NasalPair>,
    /// The left sinus, if declared.
    pub sinus_left: Option<NasalSinus>,
    /// The right sinus, if declared.
    pub sinus_right: Option<NasalSinus>,
    /// Half the distance between the nostrils on the baffle [m].
    pub nostril_x_m: f64,
    /// The nostrils' height above the mouth's centre on the baffle [m].
    pub nostril_y_m: f64,
    /// The source of the nostrils' place.
    pub nostril_source: String,
}

/// The nasal rows as read, before the card's oral sections admit them.
#[derive(Debug, Default)]
pub(super) struct NasalDraft {
    source: Option<String>,
    licence: Option<String>,
    velum: Option<(usize, f64, String)>,
    sections: Vec<TractSection>,
    pairs: Vec<NasalPair>,
    sinus: [Option<NasalSinus>; 2],
    nostrils: Option<(f64, f64, String)>,
}

fn positive(value: f64, what: &str) -> Result<f64, String> {
    if value > 0.0 { Ok(value) } else { Err(format!("{what} must be positive")) }
}

impl NasalDraft {
    /// Read one nasal row (`key` is one of [`NASAL_KEYS`]).
    pub(super) fn row(&mut self, key: &str, rest: &[&str]) -> Result<(), String> {
        match key {
            "nasal_source" => {
                arity(rest, 1, key)?;
                once(&mut self.source, rest[0].to_string(), key)
            }
            "nasal_licence" => {
                let text = rest.join(" ");
                if text.is_empty() {
                    return Err("nasal_licence needs text".into());
                }
                once(&mut self.licence, text, key)
            }
            "velum" => {
                arity(rest, 3, key)?;
                let after = usize::try_from(whole(rest, 0, "velum after_section")?)
                    .map_err(|_| "velum after_section out of range".to_string())?;
                let area = positive(field(rest, 1, "velum port area")?, "velum port area")?;
                once(&mut self.velum, (after, area, rest[2].to_string()), key)
            }
            "nasal_section" => {
                arity(rest, 2, key)?;
                if self.sections.len() >= MAX_SECTIONS {
                    return Err("too many nasal_section rows".into());
                }
                self.sections.push(TractSection {
                    area_m2: positive(field(rest, 0, "nasal_section area")?, "nasal_section area")?,
                    length_m: positive(field(rest, 1, "nasal_section length")?, "nasal_section length")?,
                });
                Ok(())
            }
            "nasal_pair" => {
                arity(rest, 3, key)?;
                if self.pairs.len() >= MAX_SECTIONS {
                    return Err("too many nasal_pair rows".into());
                }
                self.pairs.push(NasalPair {
                    left_area_m2: positive(field(rest, 0, "nasal_pair left area")?, "nasal_pair left area")?,
                    right_area_m2: positive(field(rest, 1, "nasal_pair right area")?, "nasal_pair right area")?,
                    length_m: positive(field(rest, 2, "nasal_pair length")?, "nasal_pair length")?,
                });
                Ok(())
            }
            "sinus" => {
                arity(rest, 5, key)?;
                let side = match rest[0] {
                    "left" => 0,
                    "right" => 1,
                    _ => return Err("sinus side must be `left` or `right`".into()),
                };
                let sinus = NasalSinus {
                    from_nostril_m: positive(field(rest, 1, "sinus from_nostril")?, "sinus from_nostril")?,
                    ostium_area_m2: positive(field(rest, 2, "sinus ostium area")?, "sinus ostium area")?,
                    ostium_length_m: positive(field(rest, 3, "sinus ostium length")?, "sinus ostium length")?,
                    volume_m3: positive(field(rest, 4, "sinus volume")?, "sinus volume")?,
                };
                once(&mut self.sinus[side], sinus, if side == 0 { "sinus left" } else { "sinus right" })
            }
            "nostrils" => {
                arity(rest, 3, key)?;
                let x = positive(field(rest, 0, "nostrils x")?, "nostrils x")?;
                let y = field(rest, 1, "nostrils y")?;
                once(&mut self.nostrils, (x, y, rest[2].to_string()), key)
            }
            other => Err(format!("unknown key {other:?}")),
        }
    }

    /// `None` without any nasal row; otherwise every row but `sinus` is needed, the port must open before the lips,
    /// and each sinus must open inside its passage.
    pub(super) fn admit(self, oral: &[TractSection]) -> Result<Option<NasalCard>, String> {
        let any = self.source.is_some()
            || self.licence.is_some()
            || self.velum.is_some()
            || !self.sections.is_empty()
            || !self.pairs.is_empty()
            || self.sinus.iter().any(Option::is_some)
            || self.nostrils.is_some();
        if !any {
            return Ok(None);
        }
        let need = |what: &str| format!("a nasal card needs `{what}`");
        let source = self.source.ok_or_else(|| need("nasal_source"))?;
        let licence = self.licence.ok_or_else(|| need("nasal_licence"))?;
        let (after_section, port_area_m2, attach_source) = self.velum.ok_or_else(|| need("velum"))?;
        if self.sections.is_empty() {
            return Err(need("nasal_section"));
        }
        if self.pairs.is_empty() {
            return Err(need("nasal_pair"));
        }
        let (nostril_x_m, nostril_y_m, nostril_source) = self.nostrils.ok_or_else(|| need("nostrils"))?;
        if after_section + 1 >= oral.len() {
            return Err(format!(
                "velum after_section {after_section} needs an oral section after it: the card has {}",
                oral.len()
            ));
        }
        let passage_m: f64 = self.pairs.iter().map(|p| p.length_m).sum();
        for (side, sinus) in ["left", "right"].iter().zip(&self.sinus) {
            if let Some(s) = sinus
                && !(s.from_nostril_m > CUT_SNAP_M && s.from_nostril_m < passage_m - CUT_SNAP_M)
            {
                return Err(format!(
                    "sinus {side} opens {} m before its nostril, not inside the {passage_m} m of nasal_pair rows",
                    s.from_nostril_m
                ));
            }
        }
        let [sinus_left, sinus_right] = self.sinus;
        Ok(Some(NasalCard {
            source,
            licence,
            after_section,
            port_area_m2,
            attach_source,
            sections: self.sections,
            pairs: self.pairs,
            sinus_left,
            sinus_right,
            nostril_x_m,
            nostril_y_m,
            nostril_source,
        }))
    }
}

/// The port as the graph takes it: open (the card's) or shut (no nasal edge at all: the oral tract as a two-edge chain,
/// for the gate that a shut port leaves the oral lines where they were).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Velum {
    Open,
    #[cfg_attr(not(test), allow(dead_code))]
    Shut,
}

fn cylinders(sections: impl IntoIterator<Item = (f64, f64)>) -> Vec<Segment> {
    sections
        .into_iter()
        .map(|(area, length)| Segment::Cylinder { radius: (area / core::f64::consts::PI).sqrt(), length })
        .collect()
}

/// `chain` cut at `at` metres from its start: a section the cut falls inside becomes two of the same area. A cut within
/// [`CUT_SNAP_M`] of a section end is that end, so rounding leaves no sliver cylinder.
fn split_chain(chain: &[(f64, f64)], at: f64) -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
    let (mut back, mut front) = (Vec::new(), Vec::new());
    let mut start = 0.0;
    for &(area, length) in chain {
        let end = start + length;
        if end <= at + CUT_SNAP_M {
            back.push((area, length));
        } else if start >= at - CUT_SNAP_M {
            front.push((area, length));
        } else {
            back.push((area, at - start));
            front.push((area, end - at));
        }
        start = end;
    }
    (back, front)
}

fn node(id: &str, boundary: NodeBoundary) -> NetworkNode {
    NetworkNode { id: id.into(), boundary }
}

fn port(label: &str) -> NetworkNode {
    node(label, NodeBoundary::Port { label: label.into() })
}

fn edge(id: &str, from: usize, to: usize, segments: Vec<Segment>) -> NetworkEdge {
    NetworkEdge { id: id.into(), from, to, duct: Duct { segments }, corrections: JunctionCorrections::Disabled }
}

/// The graph of the module doc and its radiating ports, with each nostril's radius.
fn nasal_graph(
    oral: &[TractSection],
    card: &NasalCard,
    termination: Termination,
    velum: Velum,
) -> (DuctNetwork, Vec<PortTermination>, [f64; 2]) {
    let radius = |area: f64| (area / core::f64::consts::PI).sqrt();
    let radiating = |label: &str, area: f64| PortTermination {
        label: label.into(),
        load: NetworkLoad::Radiation { radius: radius(area), termination },
    };
    let k = card.after_section;
    let mut nodes = vec![port("glottis"), node("velum", NodeBoundary::Junction), port("lips")];
    let mut edges = vec![
        edge("oral-back", 0, 1, cylinders(oral[..=k].iter().map(|s| (s.area_m2, s.length_m)))),
        edge("oral-front", 1, 2, cylinders(oral[k + 1..].iter().map(|s| (s.area_m2, s.length_m)))),
    ];
    let mut terminations = vec![radiating("lips", oral[oral.len() - 1].area_m2)];
    let mut nostril_radius = [0.0; 2];
    if velum == Velum::Open {
        nodes.push(node("choana", NodeBoundary::Junction));
        let choana = nodes.len() - 1;
        let passage = card
            .sections
            .iter()
            .enumerate()
            .map(|(i, s)| (if i == 0 { card.port_area_m2 } else { s.area_m2 }, s.length_m));
        edges.push(edge("nasal", 1, choana, cylinders(passage)));
        let sides = [("nostril-left", card.sinus_left), ("nostril-right", card.sinus_right)];
        for (side, (label, sinus)) in sides.into_iter().enumerate() {
            let chain: Vec<(f64, f64)> = card
                .pairs
                .iter()
                .map(|p| (if side == 0 { p.left_area_m2 } else { p.right_area_m2 }, p.length_m))
                .collect();
            let last_area = chain[chain.len() - 1].0;
            nodes.push(port(label));
            let nostril = nodes.len() - 1;
            match sinus {
                None => edges.push(edge(&format!("{label}-passage"), choana, nostril, cylinders(chain))),
                Some(s) => {
                    let length: f64 = chain.iter().map(|&(_, l)| l).sum();
                    let (back, front) = split_chain(&chain, length - s.from_nostril_m);
                    nodes.push(node(&format!("{label}-ostium"), NodeBoundary::Junction));
                    let ostium = nodes.len() - 1;
                    nodes.push(node(
                        &format!("{label}-sinus"),
                        NodeBoundary::Termination(NetworkLoad::Cavity { volume_m3: s.volume_m3 }),
                    ));
                    let cavity = nodes.len() - 1;
                    edges.push(edge(&format!("{label}-passage-back"), choana, ostium, cylinders(back)));
                    edges.push(edge(&format!("{label}-passage-front"), ostium, nostril, cylinders(front)));
                    edges.push(edge(
                        &format!("{label}-ostium"),
                        ostium,
                        cavity,
                        cylinders([(s.ostium_area_m2, s.ostium_length_m)]),
                    ));
                }
            }
            terminations.push(radiating(label, last_area));
            nostril_radius[side] = radius(last_area);
        }
    }
    (DuctNetwork { nodes, edges }, terminations, nostril_radius)
}

/// The nasal card's lines (module doc): the island's line, and the lips' and the nostrils' flow lines with their
/// receipts. With [`Velum::Shut`] the nostril lines are absent.
pub(super) struct NasalLines {
    pub(super) line: DelayedFilter,
    pub(super) inlet_radius_m: f64,
    pub(super) lips: (DelayedFilter, MouthFlowRealization),
    pub(super) nostrils: Option<[(DelayedFilter, MouthFlowRealization); 2]>,
    pub(super) nostril_radius_m: [f64; 2],
    /// The graph's regional round trip in samples, as the line's ladder seeded it (read by the passivity gate).
    #[cfg_attr(not(test), allow(dead_code))]
    pub(super) round_trip_samples: f64,
}

/// Realize [`NasalLines`] for `oral` and `card` in `gas`, every edge on `wall`, on the mechanics clock.
///
/// # Errors
/// The network owner's or the ladder's refusal, by name.
pub(super) fn nasal_lines(
    oral: &[TractSection],
    card: &NasalCard,
    gas: &GasState,
    wall: Option<&WallPin>,
    termination: Termination,
    sample_rate_hz: u32,
    velum: Velum,
) -> Result<NasalLines, String> {
    let (network, terminations, nostril_radius_m) = nasal_graph(oral, card, termination, velum);
    let gases = network.edges.iter().map(|e| (e.id.clone(), GasRecord::Legacy(*gas))).collect();
    let media = NetworkMedia::new(&network, gases, GasRecord::Legacy(*gas))
        .map_err(|e| format!("nasal network gas refused: {e:?}"))?;
    let media = match wall {
        Some(pin) => media
            .with_edge_walls(&network, network.edges.iter().map(|e| (e.id.clone(), *pin)).collect())
            .map_err(|e| format!("nasal network wall refused: {e:?}"))?,
        None => media,
    };
    let inlet_radius_m = (oral[0].area_m2 / core::f64::consts::PI).sqrt();
    let (zc, _) = glottal_port_impedance(gas, inlet_radius_m);
    let schema = crate::port_adoption::schema("voice-nasal", "glottis", PortKind::FluidPressureFlux);
    let point = NetworkDrivingPoint::new(network, "glottis".into(), terminations, schema)
        .map_err(|e| format!("nasal network refused: {e:?}"))?;
    // The island's line takes the island's ladder. The flow lines take the series mouth line's default ladder with
    // NASAL_FLOW_MAX_TAPS: the rungs and tolerances are the series line's, so a flow the series cap holds stops on the
    // same rung.
    let line_ladder = LineLadder { max_taps: GLOTTAL_LINE_MAX_TAPS, ..LineLadder::default() };
    let flow_ladder = LineLadder { max_taps: NASAL_FLOW_MAX_TAPS, ..LineLadder::default() };
    let round_trip_samples = point
        .regional_round_trip_seconds(&media)
        .map_err(|e| format!("nasal network round trip refused: {e:?}"))?
        * f64::from(sample_rate_hz);
    let (line, _) = point
        .realize_line_regional(&media, LossModel::Bessel, sample_rate_hz, zc, &line_ladder, None)
        .map_err(|e| format!("nasal tract line refused: {e}"))?;
    let flow = |label: &str| -> Result<(DelayedFilter, MouthFlowRealization), String> {
        let (ir, receipt) = point
            .realize_port_flow_ir_regional(&media, LossModel::Bessel, label, sample_rate_hz, zc, &flow_ladder)
            .map_err(|e| format!("nasal {label} flow refused: {e}"))?;
        let line = DelayedFilter::from_impulse_response(1.0 / f64::from(sample_rate_hz), ir)
            .map_err(|e| format!("nasal {label} flow line refused: {e:?}"))?;
        Ok((line, receipt))
    };
    let lips = flow("lips")?;
    let nostrils = match velum {
        Velum::Open => Some([flow("nostril-left")?, flow("nostril-right")?]),
        Velum::Shut => None,
    };
    Ok(NasalLines { line, inlet_radius_m, lips, nostrils, nostril_radius_m, round_trip_samples })
}

#[cfg(test)]
mod tests;
