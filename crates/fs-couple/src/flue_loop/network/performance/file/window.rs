//! Optional `window none|wall-cut AXIAL_CELLS` after the mouth records.
//!
//! `wall-cut` computes the mouth window's mass from its 3-D geometry
//! (`fs_bem::wall_window`) and replaces the legacy `0.61 a_eq` inertance with
//! it. The geometry is the card's own:
//! - the inlet section is the window's channel through the wall: its length is
//!   the wall thickness, and its area must equal flue width times cut-up;
//! - the window spans the flue width as arc on the bore and the cut-up along
//!   it, cut radially through the wall;
//! - the junction at the channel's inner end joins the bore and, optionally, a
//!   closed pocket: one section ending in a rigid termination (reflection 1),
//!   whose length runs from the window centre to the closed end. Pocket and
//!   bore must share one radius. Without a pocket, the bore section starts at
//!   the window centre with a flush closed end at the window's upstream edge.
//!
//! The mass is the order-4/3 extrapolation of AXIAL_CELLS and twice that,
//! with free space around the outer cylinder. The graph keeps its inlet
//! section, and the mouth inertance becomes the rest of the branch mass, so
//! the static mass from junction to far field is the computed one. The jet's
//! source coefficient and `mouth_loss_rel` are unchanged. A compact mouth
//! reaction refuses this row. Not represented: a windway block, labium or flat
//! baffle around the window, a tapered pocket, viscous edge loss.

use super::{FlueInputError, Rows, bad};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::flue_loop::{FlueError, FlueGeometry};
use fs_bem::wall_window::{
    WallWindow, WallWindowLadder, WallWindowOptions, WindowOutside, WindowUpstream,
    wall_window_ladder,
};

/// The window derived from a card and its solved mass.
#[derive(Clone, Debug, PartialEq)]
pub struct FlueWindowRecord {
    /// The solved geometry.
    pub window: WallWindow,
    /// The coarse rung's cells: axial, arc and wall layers. The fine rung
    /// doubles each.
    pub cells: [usize; 3],
    /// Both rungs and the extrapolated mass per density.
    pub ladder: WallWindowLadder,
    /// The inlet section's own `L / S` [1/m], which the graph keeps.
    pub riser_mass_per_density: f64,
}

/// Read the optional row: `Some(AXIAL_CELLS)` for `wall-cut`.
pub(super) fn read(rows: &mut Rows<'_>) -> Result<Option<usize>, FlueInputError> {
    if !rows.next_is("window") {
        return Ok(None);
    }
    let mut row = rows.row("window")?;
    let cells = match row.fields.next() {
        Some("none") => None,
        Some("wall-cut") => {
            let cells: usize = row.number()?;
            if !(2..=32).contains(&cells) {
                return Err(bad(row.line, "window wall-cut needs 2..=32 axial cells"));
            }
            Some(cells)
        }
        _ => {
            return Err(bad(
                row.line,
                "window requires none or wall-cut AXIAL_CELLS",
            ));
        }
    };
    row.end()?;
    Ok(cells)
}

/// Derive the window from the graph and the flue row, and solve its mass.
pub(super) fn solve(
    graph: &TubeNetworkSpec,
    geometry: &FlueGeometry,
    axial_cells: usize,
    line: usize,
) -> Result<FlueWindowRecord, FlueInputError> {
    let refuse = |what| Err(bad(line, what));
    let Some(inlet) = graph
        .nodes
        .iter()
        .position(|node| matches!(node, NetworkNode::Inlet))
    else {
        return refuse("window wall-cut needs the graph's inlet");
    };
    let Some(riser_index) = graph.sections.iter().position(|s| s.nodes.contains(&inlet)) else {
        return refuse("window wall-cut needs the inlet's section");
    };
    let riser = graph.sections[riser_index];
    let junction = if riser.nodes[0] == inlet {
        riser.nodes[1]
    } else {
        riser.nodes[0]
    };
    let riser_area = core::f64::consts::PI * riser.radius_m * riser.radius_m;
    if (riser_area / (geometry.flue_width_m * geometry.cut_up_m) - 1.0).abs() > 1e-3 {
        return refuse(
            "window wall-cut: the inlet section must be the window's channel, of area flue width x cut-up",
        );
    }
    if !matches!(graph.nodes[junction], NetworkNode::Junction) {
        return refuse("window wall-cut: the channel must end at a plain junction");
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
                "window wall-cut: the channel's junction must join the bore and at most one closed pocket",
            );
        }
    };
    if pocket.is_some_and(|p| p.radius_m.to_bits() != bore.radius_m.to_bits()) {
        return refuse(
            "window wall-cut: a stepped pocket is not represented; pocket and bore must share one radius",
        );
    }
    let window = WallWindow {
        bore_radius: bore.radius_m,
        wall_thickness: riser.length_m,
        arc_width: geometry.flue_width_m,
        axial_length: geometry.cut_up_m,
        upstream: WindowUpstream::ClosedEnd {
            distance: pocket.map_or(0.5 * geometry.cut_up_m, |p| p.length_m),
        },
        outside: WindowOutside::Free,
    };
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
    let options = WallWindowOptions {
        workers,
        ..WallWindowOptions::rung(&window, axial_cells)
    };
    let ladder = wall_window_ladder(&window, &options).map_err(|error| {
        FlueInputError::Physics(FlueError::Physics(AcousticRealizeError::Nonlinear(
            format!("window wall-cut: {error}"),
        )))
    })?;
    Ok(FlueWindowRecord {
        window,
        cells: [options.axial_cells, options.arc_cells, options.wall_layers],
        ladder,
        riser_mass_per_density: riser.length_m / riser_area,
    })
}
