//! Optional `window none|wall-cut AXIAL_CELLS [cheeks radial|parallel]
//! [ramp EDGE_HEIGHT_M ANGLE_DEG]|embouchure AXIAL_CELLS ACROSS_M ALONG_M
//! OPEN_M` after the mouth records.
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
//! - its outer end is covered but for a strip OPEN wide across it, centred
//!   ([`fs_bem::wall_window::WindowCover`]): the lower lip. A real lip covers
//!   the side nearer the player. In two dimensions a one-sided cover adds
//!   twice the centred one's length (by images), and that excess is not
//!   represented;
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
//! outer cylinder, and the player's face.

use super::{FlueInputError, Rows, bad};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::flue_loop::{FlueError, FlueGeometry};
use fs_bem::wall_window::{
    WallWindow, WallWindowLadder, WallWindowOptions, WindowCheeks, WindowCover, WindowOutside,
    WindowRamp, WindowUpstream, wall_window_ladder,
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
    },
    /// A flute's embouchure hole [m]: across by along, its outer end open
    /// over a centred strip.
    Embouchure { across: f64, along: f64, open: f64 },
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
            Some(WindowRow {
                axial_cells,
                shape: Shape::WallCut { cheeks, ramp },
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
            Some(WindowRow {
                axial_cells,
                shape: Shape::Embouchure {
                    across,
                    along,
                    open,
                },
            })
        }
        _ => {
            return Err(bad(
                row.line,
                "window requires none, wall-cut AXIAL_CELLS [cheeks radial|parallel] [ramp EDGE_HEIGHT_M ANGLE_DEG] or embouchure AXIAL_CELLS ACROSS_M ALONG_M OPEN_M",
            ));
        }
    };
    row.end()?;
    Ok(shape)
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
    let (width, axial_length, cheeks, ramp, cover) = match shape.shape {
        Shape::WallCut { cheeks, ramp } => {
            (geometry.flue_width_m, geometry.cut_up_m, cheeks, ramp, None)
        }
        Shape::Embouchure {
            across,
            along,
            open,
        } => (
            across,
            along,
            WindowCheeks::Parallel,
            None,
            Some(WindowCover { open_width: open }),
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
    })
}
