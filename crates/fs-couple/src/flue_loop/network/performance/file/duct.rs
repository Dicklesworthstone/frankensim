//! Geometry and gas data lowered into the existing shared tube-network owner.
//! No propagation, scattering, load or material formula is repeated here.
use super::*;
use crate::bernoulli_aperture::radiation::RadiationFlange;
use crate::flue_loop::network::performance::FlueRadiatingOutlet;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationSampling};

pub(super) struct Input {
    pub graph: TubeNetworkSpec,
    pub gases: Vec<GasState>,
    pub inner: Vec<HoleInnerLength>,
    /// The `closed-hole` nodes, each with its row's line for the clock's check.
    pub closed: Vec<(ClosedHoleLump, usize)>,
    explicit_network: bool,
}

// One choice applies to every addressed viscous/thermal pair. It changes the
// shared discrete memory law, never the gas, physical spectrum or fit budget.
fn read_loss_sampling(row: &mut Row<'_>) -> Result<RelaxationSampling, FlueInputError> {
    match row.fields.next() {
        None | Some("bilinear-midpoint") => Ok(RelaxationSampling::BilinearMidpoint),
        Some("exact-held-flow") => Ok(RelaxationSampling::ExactHeldFlow),
        Some(_) => Err(bad(row.line,
            "viscothermal sampling requires bilinear-midpoint or exact-held-flow")),
    }
}

#[cfg(test)]
mod sampling_tests;

#[cfg(test)]
mod closed_hole_tests;

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
/// chimney's length. On its own it adds neither the matching volume `t_m`
/// (their eq. 12) nor the series length `t_a` (eq. 10); the words `matching`
/// and `series` after the radius add them (below). Nothing is added at the
/// outer end, whose load is the outlet's own. The row needs a side hole: a chimney from a junction
/// of at least three sections to an open end (reflection -1), narrower than
/// the bore (`0 < d <= 1`).
///
/// `inner BORE_RADIUS_M matching` adds `t_m = b d (1 + 0.207 d^3) / 8` to the
/// chimney once: in their eq. 8, `Z_s = j Z_ch (k t_i + tan[k (t + t_m + t_r)])`,
/// it is chimney length. `series` adds their series impedance
/// `Z_a = j k Z_c t_a`, `t_a = -b d^2 / (1.78 tanh(1.84 t / b) + 0.94 + 0.540 d
/// + 0.285 d^2)` with t the row's chimney, in their circuit (Fig. 2, eqs. 3
/// and 9): `Z_a / 2` on each side of the junction and the shunt
/// `Z_h = Z_s - Z_a / 4`. Both parts go in, since keeping `t_a` in one element
/// of their matrix and not the others is not consistent (their lines 91-94).
/// `t_a` is negative and a negative inertance is not a passive load, so each
/// half is taken off the bore: each of the junction's two bore sections is
/// shortened by `(|t_a| / 2) (r_i / a)^2`, the same series mass in its own
/// area, and `-Z_a / 4` is chimney length `d^2 |t_a| / 4` (one gas across the
/// junction). The shortened bore also loses `S_a |t_a|` of volume, which the
/// source's series mass does not carry: that is this realisation's own
/// error. Its size, measured in a scratch lossless plane-wave model of the
/// Davies D whistle (HazyDune's `open_ta_size.py`, md5
/// 5b89910f1dff74c138038a708b09e191; 0.3 and 1.5 mm walls, two outer end
/// corrections): the first resonance sits 0.01 to 0.28 cent above their
/// circuit's, 1 to 18 % of the move their `t_a` itself makes, growing with
/// the number of open holes (one 0.01 cent, three 0.08 to 0.11, five 0.21
/// to 0.28). Given back as a parallel compliance `S_a |t_a| / (rho c^2)` at
/// the junction, the lost volume would leave under 0.006 cent; nothing adds
/// it, and this row needs a plain `junction` node, so a card cannot. The
/// recorder rows were not sized. `series` needs exactly two bore sections
/// at the junction, one gas on the hole and both of them, and each section
/// longer than its share.
/// `t_i` and `t_a` are their open-hole forms; `t_m` is the same for an open
/// and a closed hole (Lefebvre and Scavone 2012, Eq. 8, as fs-duct's
/// `hole_lengths` applies it). A `closed-hole` node takes none of these
/// open-hole forms: its own printed `t_m` and `t_a` come by its words (see
/// [`ClosedHoleLump`]).
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
    /// `t_m` [m] with `matching`, added to the chimney; `None` without it.
    pub matching_length_m: Option<f64>,
    /// `t_a` [m] (negative) with `series`; `None` without it.
    pub series_length_m: Option<f64>,
    /// With `series`: the junction's two bore sections and the length taken
    /// off each, `(|t_a| / 2) (r_i / a)^2` [m].
    pub series_bore: Option<[(usize, f64); 2]>,
    /// The chimney's length in the graph [m]:
    /// `t + t_i (+ t_m) (+ d^2 |t_a| / 4)`.
    pub network_length_m: f64,
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

    /// `t_m = b d (1 + 0.207 d^3) / 8` [m] (Garcia Mayen et al. eq. 12); `None`
    /// where [`Self::law`] is.
    #[must_use]
    pub fn matching_law(hole_radius_m: f64, bore_radius_m: f64) -> Option<f64> {
        Self::law(hole_radius_m, bore_radius_m).map(|_| {
            let d = hole_radius_m / bore_radius_m;
            hole_radius_m * d * (1.0 + 0.207 * d * d * d) / 8.0
        })
    }

    /// `t_a = -b d^2 / (1.78 tanh(1.84 t / b) + 0.94 + 0.540 d + 0.285 d^2)` [m]
    /// (Garcia Mayen et al. eq. 10), t the chimney height; `None` where
    /// [`Self::law`] is or unless `t > 0`.
    #[must_use]
    pub fn series_law(hole_radius_m: f64, bore_radius_m: f64, chimney_m: f64) -> Option<f64> {
        Self::law(hole_radius_m, bore_radius_m).filter(|_| chimney_m.is_finite() && chimney_m > 0.0)
            .map(|_| {
                let d = hole_radius_m / bore_radius_m;
                -hole_radius_m * d * d / (1.78 * fs_math::det::tanh(1.84 * chimney_m / hole_radius_m)
                    + 0.94 + 0.540 * d + 0.285 * d * d)
            })
    }
}

/// A closed side hole shorter than one sample step, from
/// `duct_node closed-hole HOLE_RADIUS_M HEIGHT_M` on the bore junction it
/// opens from. It stands in for the chimney's `duct_section` and its
/// `reflection 1` end, which a hole under one step cannot be: the fractional
/// owner (`TubeNetworkSpec::realize_fractional_waveguide`) refuses a section
/// shorter than one complete sample.
///
/// The hole is the network's existing shunt load, a passive series R-L-C from
/// the junction to zero pressure with no delay, holding the first two terms of
/// a closed uniform tube's input impedance,
/// `-j Z_h cot(k t) = 1 / (j w C) + j w L + O((k t)^3)`, `Z_h = rho c / S`,
/// `S = pi b^2`:
/// `C = S t / (rho c^2)`, the hole's volume as a compliance, and
/// `L = rho t / (3 S)`, its column's inertance in series (a third of the
/// column's mass: the flow falls linearly to zero at the closed end).
/// `R = 0`: the closed waveguide chimney it stands in for is lossless too.
/// b is the hole's radius and t its height on its axis from the bore wall (the
/// length a chimney row would give); rho and c are the gas of the junction's
/// sections, which must be one gas. The fs-vfit owner realises the load by the
/// bilinear (implicit-midpoint) rule, so its discrete impedance is
/// `R + s L + 1 / (s C)` at `s = j (2 / dt) tan(w dt / 2)`.
///
/// End corrections: none is added unless the row names them, as the closed
/// waveguide chimney takes none, so the two agree where the hole passes one
/// step. A closed hole has two in print (Lefebvre and Scavone 2012), and the
/// row takes them from fs-duct's `hole_lengths` (`HoleLengths::Printed`) with
/// the open hole's own words:
/// `duct_node closed-hole HOLE_RADIUS_M HEIGHT_M [matching] [series] BORE_RADIUS_M`,
/// either word or both, then a, the main bore's radius at the hole.
///
/// - `matching` adds the matching volume `t_m = b d (1 + 0.207 d^3) / 8`,
///   `d = b / a` (their Eq. 8, the same as an open hole's;
///   `HoleLengths::sealed_extra_m`) to the height: `L` and `C` are taken at
///   `t + t_m`, the first two terms of their Eq. 16, `-j Z_h cot k(t + t_m)`,
///   and the one-step admission reads `t + t_m`. A closed waveguide chimney
///   takes no words: written `t + t_m` long it is the same hole, and the
///   crossover test reads it so.
/// - `series` adds the closed series length of their Eq. 34,
///   `t_a = (-0.12 - 0.17 tanh(2.4 t / b)) b d^2`, negative
///   (`HoleLengths::series_m` with the hole closed, `sigma = 0`). Their arm
///   `Z_a = j k Z_c t_a` is a series mass only, and a negative mass is not a
///   passive load, so each of the junction's two bore sections is shortened by
///   `(|t_a| / 2) (r_i / a)^2`, the same series mass in its own area, and the
///   volume the two sections lose, `V_a = sum_i pi r_i^2 (|t_a| / 2) (r_i / a)^2`,
///   is given back to the hole: `C` holds `pi b^2 (t + t_m) + V_a`. The
///   junction then loses mass only, as the arm does. Without the volume the
///   shortening also took `V_a` off the junction, about the hole's own volume
///   on the D whistle's wide holes (17.7 against 18.1 mm^3 on the 3.3 mm hole),
///   and moved its closed notes 1.6 to 2.3 cents sharp of the arm (MB 25080).
///   `V_a` sits behind the hole column's mass `L`, not beside it: the shunt's
///   admittance differs from the arm's by `w^2 L (C^2 - C_0^2) / C_0`, `C_0`
///   without `V_a`, 4e-4 at 3 kHz on the test hole. The one-step admission
///   still reads `t (+ t_m)`. It needs exactly two sections at the junction,
///   each longer than its share.
///
/// Either word needs `0 < b < a`: the fits are printed for `b < a`, and
/// fs-duct keeps the engine's terms outside it, which add nothing to a closed
/// hole. Left out, with or without the words: the `-Z_a / 4` of the symmetric T,
/// the `k t_i` mass of Eq. 16 (the paper calls `t_i` negligible at low
/// frequency, and fs-duct leaves it out too), each under 0.001 cent on the D
/// whistle's closed holes in a scratch lossless model, and any closed-end
/// term. Without the words `t_m` is missing and `C` is low: on the D whistle's
/// six closed holes (b 2.47 to 3.31 mm, bore radius a 6.1 mm, t 0.3 mm) `t_m`
/// is 0.13 to 0.23 mm, 42 to 77 % of `t`, so `C` is 0.56 to 0.70 of the value
/// with it.
///
/// Admission: `b > 0`, `t > 0`, b no wider than any section at the junction,
/// at least two sections there, one gas on them, no `inner` side hole's
/// chimney on the same junction, and `t < c dt` (`t + t_m < c dt` with
/// `matching`) once the clock is known (checked after the duct is read): a
/// closed hole of one step or more is a `duct_section`. An open hole has no
/// lumped form here: its outlet load and its receiver live on its far node.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClosedHoleLump {
    /// The bore junction's node index.
    pub node: usize,
    /// b, the hole's radius [m].
    pub hole_radius_m: f64,
    /// t, the hole's height on its axis from the bore wall [m].
    pub height_m: f64,
    /// rho, the junction gas's density [kg/m^3].
    pub density_kg_m3: f64,
    /// c, the junction gas's sound speed [m/s].
    pub sound_speed_m_s: f64,
    /// `L = rho t / (3 pi b^2)` [Pa s^2/m^3], `t + t_m` in place of `t` with `matching`.
    pub inertance_pa_s2_m3: f64,
    /// `C = pi b^2 t / (rho c^2)` [m^3/Pa], `t + t_m` in place of `t` with
    /// `matching`; with `series` it also holds the bore volume `V_a` the
    /// shortened sections lose, `+ V_a / (rho c^2)`.
    pub compliance_m3_pa: f64,
    /// a, the main bore's radius at the hole, with either word [m]; `None` without them.
    pub bore_radius_m: Option<f64>,
    /// `t_m` [m] with `matching`; `None` without it.
    pub matching_length_m: Option<f64>,
    /// `t_a` [m] (negative) of Eq. 34 with `series`; `None` without it.
    pub series_length_m: Option<f64>,
    /// With `series`: the junction's two bore sections and the length taken
    /// off each, `(|t_a| / 2) (r_i / a)^2` [m].
    pub series_bore: Option<[(usize, f64); 2]>,
    /// With `series`: `V_a`, the volume the two shortened sections lose and
    /// `C` holds instead [m^3].
    pub series_volume_m3: Option<f64>,
}

impl ClosedHoleLump {
    /// The height `L` and `C` are taken at [m]: `t`, or `t + t_m` with `matching`.
    #[must_use]
    pub fn lumped_height_m(&self) -> f64 {
        self.matching_length_m.map_or(self.height_m, |t_m| self.height_m + t_m)
    }

    /// `(L, C)` for a closed hole of radius b and height t in a gas of density
    /// rho and sound speed c; `None` unless all four, L and C are finite and
    /// positive.
    #[must_use]
    pub fn law(hole_radius_m: f64, height_m: f64, density_kg_m3: f64, sound_speed_m_s: f64)
        -> Option<(f64, f64)>
    {
        [hole_radius_m, height_m, density_kg_m3, sound_speed_m_s].iter()
            .all(|value| value.is_finite() && *value > 0.0)
            .then(|| {
                let area = core::f64::consts::PI * hole_radius_m * hole_radius_m;
                (density_kg_m3 * height_m / (3.0 * area),
                    area * height_m / (density_kg_m3 * sound_speed_m_s * sound_speed_m_s))
            })
            .filter(|&(inertance, compliance)| inertance.is_finite() && inertance > 0.0
                && compliance.is_finite() && compliance > 0.0)
    }
}

impl Input {
    /// The clock's admission of hole lengths, once it is known: a `closed-hole`
    /// node must be shorter than one step (`t < c dt`, `t + t_m < c dt` with
    /// `matching`; one step or more is a `duct_section`), and an `inner` side
    /// hole's chimney must be at least one step long in the graph (an open hole
    /// has no lumped form here).
    pub(super) fn check_hole_steps(&self, mechanics_rate_hz: u32) -> Result<(), FlueInputError> {
        let dt = 1.0 / f64::from(mechanics_rate_hz);
        for (hole, line) in &self.closed {
            if hole.lumped_height_m() >= hole.sound_speed_m_s * dt {
                return Err(bad(*line, "duct_node closed-hole of one step or more (height, plus t_m with matching, >= c dt) is a duct_section"));
            }
        }
        for hole in &self.inner {
            if hole.network_length_m < self.gases[hole.section].sound_speed * dt {
                return Err(bad(0, "an open side hole shorter than one step in the graph has no lumped form (its outlet load and receiver live on its far node): use a finer clock"));
            }
        }
        Ok(())
    }

    pub fn read(rows: &mut Rows<'_>, default: GasState,
        sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>,
        gate: &CancelGate) -> Result<Self, FlueInputError>
    {
        let explicit_network = rows.next_is("network");
        let mut inner = Vec::new();
        let (nodes, sections, gases, closed) = if explicit_network {
            let mut row = rows.row("network")?;
            let node_count: usize = row.number()?;
            let section_count: usize = row.number()?;
            row.end()?;
            if !(2..=64).contains(&node_count) || !(1..=128).contains(&section_count) {
                return Err(bad(rows.line, "network requires 2..=64 nodes and 1..=128 sections"));
            }
            let mut nodes = Vec::with_capacity(node_count);
            let mut holes = Vec::new();
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
                    Some("closed-hole") => {
                        let hole_radius_m = row.scalar()?;
                        let height_m = row.scalar()?;
                        if !(hole_radius_m.is_finite() && hole_radius_m > 0.0
                            && height_m.is_finite() && height_m > 0.0) {
                            return Err(bad(row.line, "duct_node closed-hole needs a positive hole radius and height"));
                        }
                        // The open hole's words, then the bore radius they need ([`ClosedHoleLump`]).
                        let matching = row.fields.clone().next() == Some("matching");
                        if matching { row.fields.next(); }
                        let series = row.fields.clone().next() == Some("series");
                        if series { row.fields.next(); }
                        let line = row.line;
                        let bore_radius_m = if matching || series {
                            if row.fields.clone().next().is_none() {
                                return Err(bad(line, "duct_node closed-hole matching or series needs the main bore's radius after the words"));
                            }
                            Some(row.scalar()?)
                        } else if row.fields.clone().next().is_some() {
                            return Err(bad(line, "duct_node closed-hole takes a bore radius only after matching or series"));
                        } else { None };
                        holes.push((nodes.len(), hole_radius_m, height_m, matching, series, bore_radius_m, line));
                        // A junction until its sections and gas are read; then
                        // the hole's shunt load ([`ClosedHoleLump`]).
                        NetworkNode::Junction
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
                    let bore_radius_m = row.scalar()?;
                    let matching = row.fields.clone().next() == Some("matching");
                    if matching { row.fields.next(); }
                    let series = row.fields.clone().next() == Some("series");
                    if series { row.fields.next(); }
                    requests.push((index, bore_radius_m, matching, series, row.line));
                }
                let gas = if row.fields.clone().next() == Some("gas") {
                    row.fields.next();
                    row.gas(sources, records)?
                } else { default };
                row.end()?;
                sections.push(section);
                gases.push(gas);
            }
            // Closed holes first, so that an `inner` chimney on the same
            // junction is refused by name and never read as an open hole's.
            let mut closed = Vec::with_capacity(holes.len());
            for (node, hole_radius_m, height_m, matching, series, bore_radius_m, line) in holes {
                let incident: Vec<usize> = (0..sections.len())
                    .filter(|&j| sections[j].nodes.contains(&node)).collect();
                if incident.len() < 2 {
                    return Err(bad(line, "duct_node closed-hole needs a bore junction of at least two sections"));
                }
                if requests.iter().any(|request| incident.contains(&request.0)) {
                    return Err(bad(line, "duct_node closed-hole cannot share its junction with an inner side hole's chimney"));
                }
                if incident.iter().any(|&j| hole_radius_m > sections[j].radius_m) {
                    return Err(bad(line, "duct_node closed-hole needs a hole no wider than its junction's sections"));
                }
                let gas = gases[incident[0]];
                if incident.iter().any(|&j| gases[j] != gas) {
                    return Err(bad(line, "duct_node closed-hole needs one gas on its junction's sections"));
                }
                if bore_radius_m.is_some_and(|a| !(a.is_finite() && hole_radius_m < a)) {
                    return Err(bad(line, "duct_node closed-hole matching or series needs a hole narrower than the bore: 0 < hole radius < bore radius"));
                }
                let printed = fs_duct::hole_lengths::HoleLengths::Printed;
                let matching_length_m = bore_radius_m.filter(|_| matching)
                    .map(|a| printed.sealed_extra_m(hole_radius_m, a));
                let lumped_height_m = matching_length_m.map_or(height_m, |t_m| height_m + t_m);
                let Some((inertance, compliance)) =
                    ClosedHoleLump::law(hole_radius_m, lumped_height_m, gas.density, gas.sound_speed) else {
                    return Err(bad(line, "duct_node closed-hole needs a finite positive inertance and compliance"));
                };
                let (mut series_length_m, mut series_bore, mut series_volume_m3) = (None, None, None);
                if let Some(a) = bore_radius_m.filter(|_| series) {
                    let &[first, second] = incident.as_slice() else {
                        return Err(bad(line, "duct_node closed-hole series needs a junction of exactly two bore sections"));
                    };
                    // Eq. 34 at the physical height t: the hole closed (sigma 0), no open neighbour's mutual length.
                    let t_a = printed.series_m(hole_radius_m, height_m, a, 0.0, 0.0);
                    let mut shares = [(first, 0.0), (second, 0.0)];
                    let mut volume = 0.0;
                    for (j, share) in &mut shares {
                        let ratio = sections[*j].radius_m / a;
                        *share = -t_a / 2.0 * ratio * ratio;
                        if sections[*j].length_m <= *share {
                            return Err(bad(line, "duct_node closed-hole series needs each bore section at the hole longer than its share of t_a"));
                        }
                        sections[*j].length_m -= *share;
                        volume += core::f64::consts::PI * sections[*j].radius_m * sections[*j].radius_m * *share;
                    }
                    (series_length_m, series_bore, series_volume_m3) = (Some(t_a), Some(shares), Some(volume));
                }
                // The arm is mass only: the volume the shortened sections lose goes back into the hole's compliance.
                let compliance = series_volume_m3.map_or(compliance,
                    |volume| compliance + volume / (gas.density * gas.sound_speed * gas.sound_speed));
                let load = SeriesImpedanceSpec { resistance_pa_s_m3: 0.0, inertance_pa_s2_m3: inertance,
                    compliance_m3_pa: Some(compliance) };
                load.validate().map_err(|_| bad(line, "invalid passive closed-hole load"))?;
                let load = RelaxationImpedanceSpec::new(load, &[])
                    .map_err(|_| bad(line, "invalid passive closed-hole load"))?;
                nodes[node] = NetworkNode::Shunt { load };
                closed.push((ClosedHoleLump { node, hole_radius_m, height_m, density_kg_m3: gas.density,
                    sound_speed_m_s: gas.sound_speed, inertance_pa_s2_m3: inertance,
                    compliance_m3_pa: compliance, bore_radius_m, matching_length_m, series_length_m,
                    series_bore, series_volume_m3 }, line));
            }
            let chimneys: Vec<usize> = requests.iter().map(|request| request.0).collect();
            for (index, bore_radius_m, matching, series, line) in requests {
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
                let mut network_length_m = length_m + inner_length_m;
                let matching_length_m = matching
                    .then(|| HoleInnerLength::matching_law(radius_m, bore_radius_m)).flatten();
                if let Some(t_m) = matching_length_m { network_length_m += t_m; }
                let (mut series_length_m, mut series_bore) = (None, None);
                if series {
                    let bore: Vec<usize> = (0..sections.len())
                        .filter(|&j| j != index && sections[j].nodes.contains(&near)).collect();
                    let &[first, second] = bore.as_slice() else {
                        return Err(bad(line, "duct_section inner series needs a junction of exactly three sections: the hole and two bore sections"));
                    };
                    if chimneys.contains(&first) || chimneys.contains(&second) {
                        return Err(bad(line, "duct_section inner series needs its junction's other two sections to be bore, not another hole"));
                    }
                    if gases[first] != gases[index] || gases[second] != gases[index] {
                        return Err(bad(line, "duct_section inner series needs one gas on the hole and its two bore sections"));
                    }
                    let Some(t_a) = HoleInnerLength::series_law(radius_m, bore_radius_m, length_m) else {
                        return Err(bad(line, "duct_section inner needs a hole no wider than the bore: 0 < hole radius / bore radius <= 1"));
                    };
                    let mut shares = [(first, 0.0), (second, 0.0)];
                    for (j, share) in &mut shares {
                        let ratio = sections[*j].radius_m / bore_radius_m;
                        *share = -t_a / 2.0 * ratio * ratio;
                        if sections[*j].length_m <= *share {
                            return Err(bad(line, "duct_section inner series needs each bore section at the hole longer than its share of t_a"));
                        }
                        sections[*j].length_m -= *share;
                    }
                    let d = radius_m / bore_radius_m;
                    network_length_m += d * d * -t_a / 4.0;
                    (series_length_m, series_bore) = (Some(t_a), Some(shares));
                }
                sections[index].length_m = network_length_m;
                inner.push(HoleInnerLength { section: index, hole_radius_m: radius_m, bore_radius_m,
                    chimney_m: length_m, inner_length_m, matching_length_m, series_length_m,
                    series_bore, network_length_m });
            }
            (nodes, sections, gases, closed)
        } else {
            let mut row = rows.row("tube")?;
            let length_m = row.scalar()?;
            let radius_m = row.scalar()?;
            let max_length_error_m = row.scalar()?;
            let reflection = row.scalar()?;
            row.end()?;
            (vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
                vec![TubeSection { nodes: [0, 1], length_m, radius_m, max_length_error_m }],
                vec![default], Vec::new())
        };
        let mut graph = TubeNetworkSpec { nodes, sections,
            sound_speed_m_s: default.sound_speed, max_wave_memory_bytes: 0 };
        let inlet = graph.inlet_section_index().map_err(|source|
            FlueInputError::Physics(FlueError::Physics(source)))?;
        graph.sound_speed_m_s = gases[inlet].sound_speed;
        graph.validate_section_gases(&gases).map_err(|source|
            FlueInputError::Physics(FlueError::Physics(source)))?;
        Ok(Self { graph, gases, inner, closed, explicit_network })
    }

    pub fn read_losses(&self, rows: &mut Rows<'_>)
        -> Result<(Vec<ViscothermalSectionSpec>, RelaxationSampling), FlueInputError>
    {
        let mut row = rows.row("viscothermal")?;
        if row.fields.clone().next() == Some("none") {
            row.fields.next();
            let sampling = read_loss_sampling(&mut row)?;
            row.end()?;
            return Ok((Vec::new(), sampling));
        }
        if row.fields.clone().next() == Some("sections") {
            row.fields.next();
            let count: usize = row.number()?;
            let sampling = read_loss_sampling(&mut row)?;
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
            Ok((selections, sampling))
        } else {
            if self.explicit_network {
                return Err(bad(row.line, "network losses require none or sections COUNT with addressed loss_section records"));
            }
            let loss = ViscothermalSectionSpec { section: 0, cells: row.number()?,
                minimum_frequency_hz: row.scalar()?, maximum_frequency_hz: row.scalar()? };
            let sampling = read_loss_sampling(&mut row)?;
            row.end()?;
            Ok((vec![loss], sampling))
        }
    }

    /// The outlet rows: `outlets PRIMARY ADDITIONAL [unflanged]`, then per additional opening
    /// `outlet NODE X_M Y_M [unflanged]`. The optional trailing word selects the
    /// compact unflanged-pipe REACTION for that opening (Garcia Mayen 2021, Eq. 14's
    /// "tube without flange", via Silva et al. 2009); it needs `radiation compact-baffled`
    /// (`compact_reaction`) because only that selection installs a per-outlet load. On the
    /// `outlets` row it applies to the PRIMARY opening, a bare open end such as a whistle's
    /// bore end; on an `outlet` row, to that additional opening. The microphone model of an
    /// opening stays the baffled half-space monopole. The returned flange list is empty when no
    /// row selects one, and otherwise has one entry per opening in observed order (primary first).
    pub fn read_outlets(&self, rows: &mut Rows<'_>, compact_reaction: bool)
        -> Result<(usize, Vec<FlueRadiatingOutlet>, Vec<RadiationFlange>), FlueInputError>
    {
        if !rows.next_is("outlets") && !self.explicit_network { return Ok((1, Vec::new(), Vec::new())); }
        let mut row = rows.row("outlets")?;
        let primary: usize = row.number()?;
        let additional: usize = row.number()?;
        let primary_flange = if row.fields.clone().next() == Some("unflanged") {
            row.fields.next();
            if !compact_reaction {
                return Err(bad(row.line,
                    "outlets unflanged selects a compact open-end reaction; declare radiation compact-baffled"));
            }
            RadiationFlange::Unflanged
        } else { RadiationFlange::Baffled };
        row.end()?;
        if primary >= self.graph.nodes.len()
            || additional >= crate::pcm_wav::baffled::MAX_RAYLEIGH_INPUTS {
            return Err(bad(rows.line, "outlets exceed the physical node or receiver input budget"));
        }
        let mut outlets = Vec::with_capacity(additional);
        let mut selected = Vec::with_capacity(additional);
        for _ in 0..additional {
            let mut row = rows.row("outlet")?;
            let node: usize = row.number()?;
            if node >= self.graph.nodes.len() || node == primary
                || outlets.iter().any(|p: &FlueRadiatingOutlet| p.node == node) {
                return Err(bad(row.line, "outlet must name one distinct existing node"));
            }
            let center_m = [row.scalar()?, row.scalar()?];
            let flange = if row.fields.clone().next() == Some("unflanged") {
                row.fields.next();
                if !compact_reaction {
                    return Err(bad(row.line,
                        "outlet unflanged selects a compact open-end reaction; declare radiation compact-baffled"));
                }
                RadiationFlange::Unflanged
            } else { RadiationFlange::Baffled };
            row.end()?;
            outlets.push(FlueRadiatingOutlet { node, center_m });
            selected.push(flange);
        }
        // The shared DuctPressure owner validates terminal degree, overlap,
        // common-pressure exterior and the geometry/band of every aperture.
        let flanges = if primary_flange == RadiationFlange::Unflanged
            || selected.contains(&RadiationFlange::Unflanged) {
            core::iter::once(primary_flange).chain(selected).collect()
        } else { Vec::new() };
        Ok((primary, outlets, flanges))
    }
}
