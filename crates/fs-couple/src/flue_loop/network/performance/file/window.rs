//! Optional `window none|wall-cut AXIAL_CELLS [cheeks radial|parallel]
//! [ramp EDGE_HEIGHT_M ANGLE_DEG] [mass ernoult-fabre-2017]|embouchure
//! AXIAL_CELLS ACROSS_M ALONG_M OPEN_M [lip centred|one-sided [face SETBACK_M]]`
//! after the mouth records, then optional `jet_forcing_area open-strip [source]`
//! on an embouchure.
//!
//! `mass ernoult-fabre-2017` takes the window's mass from Ernoult and Fabre,
//! Acta Acustica united with Acustica 103 (2017) 106-116, in place of the
//! computed one (see [`ErnoultFabreWindow`]). Their closed form was fitted to
//! 444 finite-element windows (standard deviation under 5 %) and their FE
//! agrees with measured recorder windows to about 10 %. The computed solve
//! still runs, and its value is kept in the record for comparison. Only its
//! inner end correction has been checked against a cited law; its outer side
//! is untested. The law needs parallel cheeks and a closed end flush with the
//! window's upstream edge. Outside the spans of their FE set it refuses.
//!
//! `wall-cut` computes the mouth window's mass from its 3-D geometry
//! (`fs_bem::wall_window`) and replaces the legacy `0.61 a_eq` inertance with
//! it. The geometry is the card's own:
//! - the inlet section is the window's channel through the wall: its length is
//!   the wall thickness, and its area must equal flue width times cut-up;
//! - the window is the cut-up long. With radial cheeks (the default) it spans
//!   the flue width as arc on the bore. With `cheeks parallel` its side walls
//!   stand the flue width apart, as the windway's do;
//! - `ramp` makes the downstream face a labium ramp: straight up to the edge
//!   height above the bore, then receding downstream at the angle (degrees)
//!   to the axis until it meets the outer surface. Without it the face rises
//!   straight through the wall;
//! - the junction at the channel's inner end joins the bore and, optionally, a
//!   closed pocket: one section ending in a rigid termination (reflection 1),
//!   whose length runs from the window centre to the closed end. Without a
//!   pocket, the bore section starts at the window centre with a flush closed
//!   end at the window's upstream edge.
//!
//! `embouchure` is a flute's: the window is its embouchure hole, not the jet's
//! flue width by cut-up, and the lip covers part of it.
//! - the inlet section is the hole's riser, as long as the wall the solve
//!   cuts; its area must equal ACROSS times ALONG;
//! - the hole is a straight rectangular prism, ACROSS between parallel walls
//!   (the jet's direction, around the bore) and ALONG the axis;
//! - its outer end is covered but for a strip OPEN wide across it
//!   ([`fs_bem::wall_window::WindowCover`]): the lower lip. `lip centred`, the
//!   default, centres the strip. `lip one-sided` puts it against the hole's
//!   far edge, the lip covering the rest from the player's side, as a real
//!   lip does. `face SETBACK_M` after it stands the player's face in the free
//!   space outside: a rigid half-plane through the bore's axis, rising from
//!   the outer surface SETBACK behind the lip's edge (as arc on that
//!   surface) along the whole bore;
//! - pocket and bore are the cork cavity and the head joint, as above.
//!
//! Pocket and bore must be one bore: radii within `POCKET_RADIUS_TOLERANCE`
//! (2.5 %) of each other, solved at the bore section's radius. A card gives
//! each section the radius of the cylinder of equal volume, so a tapering bore
//! gives the short pocket and the bore section beside it radii that differ by
//! the taper between their volume centres, without a step in the instrument:
//! 1.9 % for the flute's cork cavity and first head-joint section. The network
//! keeps both sections and radii; only the window's near field sees one bore.
//!
//! A radial prism's mass is the order-4/3 extrapolation of AXIAL_CELLS and
//! twice that; any other cut is a staircase, extrapolated from AXIAL_CELLS,
//! twice and four times that at the order the three show. Free space lies
//! around the outer cylinder. The computed mass holds the channel, so the
//! graph keeps the inlet section only as the inlet line it needs: one transit
//! sample of its gas at the mechanics clock, the shortest the fractional image
//! admits. A line as long as the channel rings against the window mass, yet a
//! channel as wide as the window carries plane waves only below its first
//! cross mode, `c / 2W` (13.5 kHz across the recorder's window); the
//! recorder's 6.575 mm wall rang so at 18 to 21 kHz. The mouth
//! inertance becomes the rest of the branch mass, so the static mass from
//! junction to far field is the computed one. The jet's source coefficient and
//! `mouth_loss_rel` are unchanged. A compact mouth reaction refuses this row.
//! Not represented: the channel air's compressibility, a windway block or flat
//! baffle around the window, a tapered pocket, viscous edge loss; for an
//! embouchure, the hole's taper and undercut, a lip plate other than the
//! outer cylinder, the lip's thickness and shape, and a face other than a
//! half-plane through the axis (inclined, curved, of finite extent, or with
//! the lips' opening in it).

use super::{FlueInputError, Rows, bad};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::flue_loop::{FlueError, FlueGeometry};
use fs_bem::wall_window::{
    CoverSide, WallWindow, WallWindowLadder, WallWindowOptions, WindowCheeks, WindowCover,
    WindowOutside, WindowRamp, WindowUpstream, wall_window_ladder,
};
use fs_material::gas::GasState;

/// Largest relative difference between pocket and bore radii that the
/// window solve takes as one bore (see the module notes).
const POCKET_RADIUS_TOLERANCE: f64 = 0.025;

/// A `window wall-cut` or `window embouchure` row: the coarse rung and the
/// cut's shape.
#[derive(Clone, Copy, Debug)]
pub(super) struct WindowRow {
    axial_cells: usize,
    shape: Shape,
}

#[derive(Clone, Copy, Debug)]
enum Shape {
    /// The flue's window, flue width by cut-up.
    WallCut {
        cheeks: WindowCheeks,
        ramp: Option<WindowRamp>,
        ernoult_fabre: bool,
    },
    /// A flute's embouchure hole [m]: across by along, its outer end open
    /// over a strip on the lip's `side`.
    Embouchure {
        across: f64,
        along: f64,
        open: f64,
        side: CoverSide,
    },
}

/// The window derived from a card and its solved mass.
#[derive(Clone, Debug, PartialEq)]
pub struct FlueWindowRecord {
    /// The solved geometry.
    pub window: WallWindow,
    /// The coarse rung's cells: axial, arc and wall layers. Each further rung
    /// doubles each.
    pub cells: [usize; 3],
    /// The rungs and the extrapolated mass per density.
    pub ladder: WallWindowLadder,
    /// The inlet line the graph keeps in place of the channel [m]: one
    /// transit sample at the mechanics clock.
    pub channel_line_m: f64,
    /// That line's own `L / S` [1/m].
    pub riser_mass_per_density: f64,
    /// With `mass ernoult-fabre-2017`: their closed form, which then gives
    /// the mass. `ladder` keeps the computed value for comparison.
    pub ernoult_fabre: Option<ErnoultFabreWindow>,
}

impl FlueWindowRecord {
    /// The static mass per density the network uses, from the junction at
    /// the window centre to the far field [1/m]: Ernoult and Fabre's, moved
    /// to that junction, when the row selects it; otherwise the computed one.
    #[must_use]
    pub fn mass_per_density(&self) -> f64 {
        self.ernoult_fabre
            .map_or(self.ladder.mass_per_density, |ef| ef.junction_mass_per_density)
    }
}

/// The spans of Ernoult and Fabre's 444 finite-element windows, read from
/// their paper: it gives no table of the set. Angle: section 4.2 (0.09 rad
/// is their lowest, pi/2 the right angle). `W tan(alpha) / l_e`: the abscissa
/// of their Fig. 8. `l_e / H`: Fig. 9. W: Fig. 11. H and the head radius:
/// the measured recorders of their Table 1, whose geometry the set follows.
pub const ERNOULT_FABRE_SPANS: ErnoultFabreSpans = ErnoultFabreSpans {
    angle_rad: [0.09, core::f64::consts::FRAC_PI_2],
    width_tan_over_ears: [0.015, 4.0],
    ears_over_height: [0.0, 2.5],
    width_m: [0.002, 0.010],
    height_m: [0.0075, 0.0192],
    bore_radius_m: [0.0053, 0.016],
};

/// Closed parameter intervals of a fitted law (see [`ERNOULT_FABRE_SPANS`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ErnoultFabreSpans {
    /// Ramp angle to the axis [rad].
    pub angle_rad: [f64; 2],
    /// `W tan(alpha) / l_e`, checked only below the right angle.
    pub width_tan_over_ears: [f64; 2],
    /// `l_e / H`.
    pub ears_over_height: [f64; 2],
    /// Window length W, channel exit to edge [m].
    pub width_m: [f64; 2],
    /// Window width H, between the cheeks [m].
    pub height_m: [f64; 2],
    /// Head bore radius r_p [m].
    pub bore_radius_m: [f64; 2],
}

/// Ernoult and Fabre 2017's window inertance, per density [1/m], with its
/// reference plane. Their eq. 18: `M_w = m_ra + m_alpha + rho l_r / (W H)`,
/// at the EDGE plane: their head length runs "between the outlet of the head
/// and the edge". `m_ra` is their eq. 16, the right-angle bend on a
/// rectangular tube of the bore's area (`delta = 2 W H / S_p`). `m_alpha` is
/// their eq. 23, the ears rising at alpha from the edge, or eq. 18's
/// `rho l_e / (W H)` at the right angle (no `ramp`). `l_r = 0.695 sqrt(W H / pi)`.
/// The card's junction sits at the window centre, W/2 upstream of the edge,
/// and its bore line starts there, so that W/2 of bore is already in the
/// network: the junction's mass is the edge's minus `(W/2) / S_p`.
/// Their model has the edge in the window plane; a `ramp`'s straight rise
/// below the edge (`edge_height_m`) is not in it, and is not added.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ErnoultFabreWindow {
    /// W: window length, channel exit to edge (the cut-up) [m].
    pub width_m: f64,
    /// H: window width between the parallel cheeks (the flue width) [m].
    pub height_m: f64,
    /// l_e: height of the ears, the wall the window cuts [m].
    pub ears_m: f64,
    /// alpha: the ramp's angle to the axis, pi/2 without a ramp [rad].
    pub angle_rad: f64,
    /// r_p: the head bore's radius [m].
    pub bore_radius_m: f64,
    /// `m_ra / rho`, eq. 16 [1/m].
    pub right_angle_per_density: f64,
    /// `m_alpha / rho`, eq. 23 (eq. 18 at the right angle) [1/m].
    pub ears_per_density: f64,
    /// `l_r`, the radiation length [m].
    pub radiation_length_m: f64,
    /// `M_w / rho` at the edge plane [1/m].
    pub edge_mass_per_density: f64,
    /// `(W/2) / S_p`: the bore from the window centre to the edge [1/m].
    pub centre_to_edge_per_density: f64,
    /// `M_w / rho` at the card's junction, the window centre [1/m].
    pub junction_mass_per_density: f64,
}

impl ErnoultFabreWindow {
    /// Evaluate the law for one window, without the span check.
    /// `angle_rad` is `None` for the right angle (no ramp).
    #[must_use]
    pub fn new(width_m: f64, height_m: f64, ears_m: f64, angle_rad: Option<f64>,
        bore_radius_m: f64) -> Self
    {
        use core::f64::consts::PI;
        let (w, h, le, rp) = (width_m, height_m, ears_m, bore_radius_m);
        let bore_area = PI * rp * rp;
        let window_area = w * h;
        let delta = 2.0 * window_area / bore_area;
        let right_angle_per_density = 2.0 / (PI * h)
            * (fs_math::det::ln((1.0 + 0.25 * delta * delta) / (2.0 * delta))
                + (2.0 / delta) * (1.0 - 0.25 * delta * delta) * fs_math::det::atan(0.5 * delta));
        let ears_per_density = match angle_rad {
            None => le / window_area,
            Some(alpha) => {
                let t = fs_math::det::tan(alpha);
                (t * fs_math::det::ln(1.0 + le / (w * t))
                    + fs_math::det::exp(-alpha * fs_math::det::ln(32.0)) * 0.85
                        * fs_math::det::sqrt(le / h)) / h
            }
        };
        let radiation_length_m = 0.695 * fs_math::det::sqrt(window_area / PI);
        let edge_mass_per_density =
            right_angle_per_density + ears_per_density + radiation_length_m / window_area;
        let centre_to_edge_per_density = 0.5 * w / bore_area;
        Self {
            width_m, height_m, ears_m, angle_rad: angle_rad.unwrap_or(core::f64::consts::FRAC_PI_2),
            bore_radius_m, right_angle_per_density, ears_per_density, radiation_length_m,
            edge_mass_per_density, centre_to_edge_per_density,
            junction_mass_per_density: edge_mass_per_density - centre_to_edge_per_density,
        }
    }

    /// The refusal for the first parameter outside [`ERNOULT_FABRE_SPANS`], if any.
    #[must_use]
    pub fn outside_spans(&self) -> Option<&'static str> {
        let s = ERNOULT_FABRE_SPANS;
        let inside = |v: f64, [lo, hi]: [f64; 2]| v >= lo && v <= hi;
        let right = self.angle_rad.to_bits() == core::f64::consts::FRAC_PI_2.to_bits();
        if !inside(self.angle_rad, s.angle_rad) {
            return Some("window mass ernoult-fabre-2017: the ramp angle lies outside their finite-element set");
        }
        if !right && !inside(self.width_m * fs_math::det::tan(self.angle_rad) / self.ears_m,
            s.width_tan_over_ears) {
            return Some("window mass ernoult-fabre-2017: W tan(alpha) / l_e lies outside their finite-element set");
        }
        if !inside(self.ears_m / self.height_m, s.ears_over_height) {
            return Some("window mass ernoult-fabre-2017: l_e / H lies outside their finite-element set");
        }
        if !inside(self.width_m, s.width_m) {
            return Some("window mass ernoult-fabre-2017: the cut-up lies outside their finite-element set");
        }
        if !inside(self.height_m, s.height_m) {
            return Some("window mass ernoult-fabre-2017: the flue width lies outside their finite-element set");
        }
        if !inside(self.bore_radius_m, s.bore_radius_m) {
            return Some("window mass ernoult-fabre-2017: the head bore radius lies outside their finite-element set");
        }
        None
    }
}

/// Read the optional row: `Some` for `wall-cut`.
pub(super) fn read(rows: &mut Rows<'_>) -> Result<Option<WindowRow>, FlueInputError> {
    if !rows.next_is("window") {
        return Ok(None);
    }
    let mut row = rows.row("window")?;
    let shape = match row.fields.next() {
        Some("none") => None,
        Some("wall-cut") => {
            let axial_cells: usize = row.number()?;
            if !(2..=32).contains(&axial_cells) {
                return Err(bad(row.line, "window wall-cut needs 2..=32 axial cells"));
            }
            let cheeks = if row.fields.clone().next() == Some("cheeks") {
                row.fields.next();
                match row.fields.next() {
                    Some("radial") => WindowCheeks::Radial,
                    Some("parallel") => WindowCheeks::Parallel,
                    _ => return Err(bad(row.line, "window cheeks requires radial or parallel")),
                }
            } else {
                WindowCheeks::Radial
            };
            let ramp = if row.fields.clone().next() == Some("ramp") {
                row.fields.next();
                let edge_height = row.scalar()?;
                let degrees = row.scalar()?;
                if !(edge_height >= 0.0 && degrees > 0.0 && degrees < 90.0) {
                    return Err(bad(
                        row.line,
                        "window ramp needs an edge height of at least 0 m and an angle inside (0, 90) degrees",
                    ));
                }
                Some(WindowRamp {
                    edge_height,
                    angle: degrees.to_radians(),
                })
            } else {
                None
            };
            let ernoult_fabre = if row.fields.clone().next() == Some("mass") {
                row.fields.next();
                match row.fields.next() {
                    Some("ernoult-fabre-2017") => true,
                    _ => return Err(bad(row.line, "window mass requires ernoult-fabre-2017")),
                }
            } else {
                false
            };
            Some(WindowRow {
                axial_cells,
                shape: Shape::WallCut { cheeks, ramp, ernoult_fabre },
            })
        }
        Some("embouchure") => {
            let axial_cells: usize = row.number()?;
            if !(2..=32).contains(&axial_cells) {
                return Err(bad(row.line, "window embouchure needs 2..=32 axial cells"));
            }
            let (across, along, open) = (row.scalar()?, row.scalar()?, row.scalar()?);
            if !(across > 0.0 && along > 0.0 && open > 0.0 && open < across) {
                return Err(bad(
                    row.line,
                    "window embouchure needs a positive hole and an open strip narrower than the hole",
                ));
            }
            let side = if row.fields.clone().next() == Some("lip") {
                row.fields.next();
                match row.fields.next() {
                    Some("centred") => CoverSide::Centred,
                    Some("one-sided") => {
                        let face = if row.fields.clone().next() == Some("face") {
                            row.fields.next();
                            let setback = row.scalar()?;
                            if !(setback >= 0.0) {
                                return Err(bad(
                                    row.line,
                                    "window face needs a setback of at least 0 m",
                                ));
                            }
                            Some(setback)
                        } else {
                            None
                        };
                        CoverSide::OneSided { face }
                    }
                    _ => return Err(bad(row.line, "window lip requires centred or one-sided")),
                }
            } else {
                CoverSide::Centred
            };
            Some(WindowRow {
                axial_cells,
                shape: Shape::Embouchure {
                    across,
                    along,
                    open,
                    side,
                },
            })
        }
        _ => {
            return Err(bad(
                row.line,
                "window requires none, wall-cut AXIAL_CELLS [cheeks radial|parallel] [ramp EDGE_HEIGHT_M ANGLE_DEG] [mass ernoult-fabre-2017] or embouchure AXIAL_CELLS ACROSS_M ALONG_M OPEN_M [lip centred|one-sided [face SETBACK_M]]",
            ));
        }
    };
    row.end()?;
    Ok(shape)
}

/// The `jet_forcing_area open-strip [source]` row.
#[derive(Clone, Copy, Debug)]
pub(super) struct JetForcingRow {
    /// The embouchure's open strip, ALONG times OPEN [m^2].
    pub(super) area_m2: f64,
    /// Whether the jet's source coefficient moves to the strip with its transverse velocity.
    pub(super) source: bool,
}

/// Read the optional `jet_forcing_area open-strip [source]` right after the window row: the
/// embouchure's open strip for the jet's transverse velocity and, with `source`, its source
/// coefficient.
pub(super) fn read_jet_forcing_area(
    rows: &mut Rows<'_>,
    window: Option<WindowRow>,
) -> Result<Option<JetForcingRow>, FlueInputError> {
    if !rows.next_is("jet_forcing_area") {
        return Ok(None);
    }
    let mut row = rows.row("jet_forcing_area")?;
    if row.fields.next() != Some("open-strip") {
        return Err(bad(row.line, "jet_forcing_area requires open-strip"));
    }
    let Some(WindowRow { shape: Shape::Embouchure { along, open, .. }, .. }) = window else {
        return Err(bad(row.line, "jet_forcing_area open-strip requires a window embouchure row before it"));
    };
    let source = match row.fields.next() {
        None => false,
        Some("source") => true,
        Some(_) => return Err(bad(row.line, "jet_forcing_area open-strip takes only the term source")),
    };
    row.end()?;
    Ok(Some(JetForcingRow { area_m2: along * open, source }))
}

/// Derive the window from the graph and the flue row, solve its mass, and
/// shorten the graph's inlet section to one transit sample.
pub(super) fn solve(
    graph: &mut TubeNetworkSpec,
    gases: &[GasState],
    mechanics_rate_hz: u32,
    geometry: &FlueGeometry,
    shape: WindowRow,
    line: usize,
) -> Result<FlueWindowRecord, FlueInputError> {
    let refuse = |what| Err(bad(line, what));
    let Some(inlet) = graph
        .nodes
        .iter()
        .position(|node| matches!(node, NetworkNode::Inlet))
    else {
        return refuse("window needs the graph's inlet");
    };
    let Some(riser_index) = graph.sections.iter().position(|s| s.nodes.contains(&inlet)) else {
        return refuse("window needs the inlet's section");
    };
    let riser = graph.sections[riser_index];
    let junction = if riser.nodes[0] == inlet {
        riser.nodes[1]
    } else {
        riser.nodes[0]
    };
    let riser_area = core::f64::consts::PI * riser.radius_m * riser.radius_m;
    let (window_area, channel) = match shape.shape {
        Shape::WallCut { .. } => (
            geometry.flue_width_m * geometry.cut_up_m,
            "window wall-cut: the inlet section must be the window's channel, of area flue width x cut-up",
        ),
        Shape::Embouchure { across, along, .. } => (
            across * along,
            "window embouchure: the inlet section must be the hole's riser, of area across x along",
        ),
    };
    if (riser_area / window_area - 1.0).abs() > 1e-3 {
        return refuse(channel);
    }
    if !matches!(graph.nodes[junction], NetworkNode::Junction) {
        return refuse("window: the channel must end at a plain junction");
    }
    let closed = |s: &&crate::bernoulli_aperture::network::TubeSection| {
        let far = if s.nodes[0] == junction {
            s.nodes[1]
        } else {
            s.nodes[0]
        };
        matches!(graph.nodes[far],
            NetworkNode::Termination { reflection } if reflection.to_bits() == 1f64.to_bits())
    };
    let others: Vec<_> = graph
        .sections
        .iter()
        .enumerate()
        .filter(|(i, s)| *i != riser_index && s.nodes.contains(&junction))
        .map(|(_, s)| s)
        .collect();
    let (bore, pocket) = match others.as_slice() {
        [bore] if !closed(bore) => (*bore, None),
        [a, b] if closed(a) != closed(b) => {
            if closed(a) {
                (*b, Some(*a))
            } else {
                (*a, Some(*b))
            }
        }
        _ => {
            return refuse(
                "window: the channel's junction must join the bore and at most one closed pocket",
            );
        }
    };
    if pocket
        .is_some_and(|p| !((p.radius_m / bore.radius_m - 1.0).abs() <= POCKET_RADIUS_TOLERANCE))
    {
        return refuse(
            "window: a stepped pocket is not represented; pocket and bore radii must agree within 2.5 %",
        );
    }
    let a = bore.radius_m;
    let ernoult_fabre = match shape.shape {
        Shape::WallCut { cheeks, ramp, ernoult_fabre: true } => {
            if cheeks != WindowCheeks::Parallel {
                return refuse("window mass ernoult-fabre-2017 needs parallel cheeks: their window is a W x H rectangle");
            }
            // Their closed end is flush with the window's upstream edge; a
            // card without a pocket has that end by construction.
            if pocket.is_some_and(|p| !((p.length_m / (0.5 * geometry.cut_up_m) - 1.0).abs() <= 1e-3)) {
                return refuse("window mass ernoult-fabre-2017 needs the closed end flush with the window's upstream edge: a pocket half the cut-up long");
            }
            let law = ErnoultFabreWindow::new(geometry.cut_up_m, geometry.flue_width_m,
                riser.length_m, ramp.map(|r| r.angle), a);
            if let Some(what) = law.outside_spans() {
                return refuse(what);
            }
            Some(law)
        }
        _ => None,
    };
    let (width, axial_length, cheeks, ramp, cover) = match shape.shape {
        Shape::WallCut { cheeks, ramp, .. } => {
            (geometry.flue_width_m, geometry.cut_up_m, cheeks, ramp, None)
        }
        Shape::Embouchure {
            across,
            along,
            open,
            side,
        } => (
            across,
            along,
            WindowCheeks::Parallel,
            None,
            Some(WindowCover {
                open_width: open,
                side,
            }),
        ),
    };
    let arc_width = match cheeks {
        WindowCheeks::Radial => width,
        WindowCheeks::Parallel => {
            if !(width < 2.0 * a) {
                return refuse(
                    "window: parallel cheeks must stand closer than the bore's diameter",
                );
            }
            2.0 * a * fs_math::det::asin(0.5 * width / a)
        }
    };
    let window = WallWindow {
        bore_radius: a,
        wall_thickness: riser.length_m,
        arc_width,
        axial_length,
        cheeks,
        ramp,
        upstream: WindowUpstream::ClosedEnd {
            distance: pocket.map_or(0.5 * axial_length, |p| p.length_m),
        },
        outside: WindowOutside::Free,
        cover,
    };
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
    let options = WallWindowOptions {
        workers,
        ..WallWindowOptions::rung(&window, shape.axial_cells)
    };
    // The same product the fractional image forms, so the line is exactly one
    // sample. A channel shorter than that never had a line.
    let channel_line_m = gases[riser_index].sound_speed * (1.0 / f64::from(mechanics_rate_hz));
    if !(channel_line_m.is_finite() && channel_line_m > 0.0 && riser.length_m >= channel_line_m) {
        return refuse("window: the channel must span at least one transit sample");
    }
    let ladder = wall_window_ladder(&window, &options).map_err(|error| {
        FlueInputError::Physics(FlueError::Physics(AcousticRealizeError::Nonlinear(
            format!("window: {error}"),
        )))
    })?;
    graph.sections[riser_index].length_m = channel_line_m;
    Ok(FlueWindowRecord {
        window,
        cells: [options.axial_cells, options.arc_cells, options.wall_layers],
        ladder,
        channel_line_m,
        riser_mass_per_density: channel_line_m / riser_area,
        ernoult_fabre,
    })
}
