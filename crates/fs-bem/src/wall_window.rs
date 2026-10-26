//! Potential-flow acoustic mass of a window cut through the wall of a
//! cylindrical bore: a flue instrument's mouth as the air column sees it.
//!
//! The bore has radius `a` and a rigid wall of thickness `t`. The window
//! meets the bore over the arc `|theta| <= theta_w` and the axial interval
//! `|x| <= L/2`, and its walls are rigid. Its side walls (cheeks) run along
//! radii or parallel to the window's mid-plane. Its downstream face either
//! rises straight through the wall or, as a labium ramp, rises straight to an
//! edge height and then recedes downstream at a fixed angle to the axis.
//! Upstream, the bore ends in a flat closed end (a flue's pocket) or in a
//! port; downstream it continues to a port. Outside, the wall is an infinite
//! rigid cylinder of radius `a + t` in free space. Alternatively the window's
//! outer face is a pressure release, which leaves no outer region and is used
//! to check the inner end correction on its own.
//!
//! The air is incompressible and irrotational. The ports carry plane flow at
//! potential 1, the far exterior sits at 0, and `Q` is the volume flux
//! through the window per unit potential difference. Extrapolate the bore's
//! plane-flow potential to the window centre:
//! - one downstream port at distance `X`, bore area `S`: `C = 1 - Q X / S`;
//! - ports at `X_d` and `X_u` on both sides: `C = 1 - Q / (S (1/X_d + 1/X_u))`.
//!
//! The window branch runs from the 1-D network's junction (on the axis, at
//! the window centre) to the far field. Its acoustic mass is `rho C / Q`,
//! which covers the channel itself, the inner (bore-side) and outer
//! corrections, and the pocket's near-field effect. Compliance stays with the
//! network, whose pocket and bore remain 1-D sections.
//!
//! Discretization: finite volumes on one cylindrical `(r, theta, x)` grid.
//! It is uniform in `theta` with `N` cells around the bore, has `arc_cells`
//! across the window and `axial_cells` along it, and grades geometrically
//! away from the window in `r` and `x`.
//! - Bore and exterior are `theta`-invariant away from the window, so each is
//!   solved exactly per discrete `theta` mode. Each mode is a
//!   block-tridiagonal system in `x` with dense radial blocks, eliminated by
//!   block Cholesky. The result is the potential on the region's window face
//!   for unit fluxes through each face cell.
//! - A prism (radial cheeks, no ramp) has a channel of `wall_layers` radial
//!   layers that is diagonal in the window cells' DCT-II basis. Per mode, it
//!   maps its two face potentials to their fluxes. The three regions couple
//!   through the two faces in that basis. Geometry and forcing are symmetric
//!   about `theta = 0`, so only the `theta`-even DCT modes enter, and each
//!   bore mode's contribution has rank one in `theta`. The coupled system has
//!   `2 ceil(arc_cells / 2) axial_cells` unknowns (half that with a release
//!   face) and is solved by LU.
//! - Any other cut is a staircase of cut cells. Each wall layer is open
//!   where the cut crosses its mid-radius: within the cheeks' angle there,
//!   and from the window's upstream edge to its downstream edge or, above the
//!   labium edge, to the ramp. A cell the cut crosses is open over that
//!   fraction of its arc or axial extent, and a cut within `1e-9` of a cell
//!   edge lies on it. Its links scale with their open face area, and a link
//!   into it ends at its open part's centroid. The bore face is the window's
//!   cells. The outer face is the last layer's, each cell open over the
//!   smaller of its own fraction and the cut's at the outer surface. Each
//!   face cell links to its channel cell over its open fraction, and each
//!   face has its own DCT bases. The exterior's grid ends every open part of
//!   the outer face on a cell edge, so no flux spreads over a closed part of
//!   a cell: axially the face's last cell ends where the opening does; in
//!   `theta` the exterior has its own grid, the bore's with the face's edge
//!   cell split, whose modes solve `K v = mu M v` with `M` the cell widths
//!   and `K` the links of one over the centre spacing. Each face cell's
//!   half-cell drop covers its open part. The face potentials are unknowns. The
//!   bore and the exterior enter as their face maps inverted, which are
//!   their Dirichlet-to-Neumann maps. The channel's cells are folded onto
//!   the `theta >= 0` half and solved with the face potentials by
//!   Jacobi-preconditioned conjugate gradients to a declared relative
//!   residual.
//! - The exterior's far boundary is a monopole Robin condition about the
//!   origin. In the cases checked, its error falls about eightfold for each
//!   doubling of the far radius.
//!
//! The arc is resolved as `arc_cells` whole cells of `2 pi / N`, with `N` the
//! nearest integer. The resolved arc can therefore differ from the requested
//! one; it is reported, and a ladder should use rungs that resolve it alike.
//! A prism's edges are re-entrant right angles, where the velocity is
//! singular like `d^(-1/3)`, so its mass converges at order 4/3 in the window
//! cell size (1.30 observed). A staircase's cut cells follow the cut
//! without rounding it to whole cells, but its edges are not all right
//! angles and its cut cells are first-order, so its order is observed from
//! three rungs, not assumed. A caller ladders the cells and extrapolates.
//! Nothing is fitted.
//!
//! Not represented: cheeks that are neither radial nor parallel; an upstream
//! face that is not straight; a flat baffle, windway or labium outside the
//! wall; a pocket that is not flat; viscosity; compressibility.

use core::f64::consts::PI;
use std::fmt;
use std::ops::Range;

use fs_la::factor::{FactorError, cholesky, lu};
use fs_math::det;

/// What closes the bore upstream of the window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WindowUpstream {
    /// A flat rigid end, as in a flue's closed pocket.
    ClosedEnd {
        /// Window centre to the end [m], at least half the window length.
        distance: f64,
    },
    /// A plane-flow port at the downstream port's potential, so that the
    /// window is a side branch of a through bore.
    Port {
        /// Window centre to the port [m], beyond the window.
        distance: f64,
    },
}

/// What lies beyond the window's outer face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowOutside {
    /// Free space around the bore's outer cylinder.
    Free,
    /// A pressure release on the outer face, with no outer correction.
    Release,
}

/// How the window's side walls run through the wall.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowCheeks {
    /// Along radii: the window spans the same angle at every radius.
    Radial,
    /// Parallel to the window's mid-plane, as a straight cutter or a
    /// moulded slot leaves them. They meet the bore at the arc `arc_width`,
    /// so they stand `2 a sin(arc_width / 2a)` apart.
    Parallel,
}

/// A labium ramp: the window's downstream face rises straight from the bore
/// to `edge_height`, then recedes downstream at `angle` to the bore axis
/// until it meets the outer surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowRamp {
    /// Height of the labium edge above the bore surface [m], in `[0, t]`.
    /// At `t` the ramp has no extent and the face rises straight through.
    pub edge_height: f64,
    /// Angle between the ramp and the bore axis [rad], in `(0, pi/2)`.
    pub angle: f64,
}

/// A window cut through the wall of a cylindrical bore.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallWindow {
    /// Bore radius `a` [m].
    pub bore_radius: f64,
    /// Wall thickness `t` at the window [m].
    pub wall_thickness: f64,
    /// Window width as arc length on the bore surface, `2 theta_w a` [m].
    pub arc_width: f64,
    /// Window length along the bore axis at the bore surface [m].
    pub axial_length: f64,
    /// Side walls.
    pub cheeks: WindowCheeks,
    /// Downstream face; `None` rises straight through the wall.
    pub ramp: Option<WindowRamp>,
    /// Upstream closure.
    pub upstream: WindowUpstream,
    /// Outer side.
    pub outside: WindowOutside,
}

impl WallWindow {
    /// Radial cheeks and no ramp: the cut whose channel is diagonal in the
    /// window's DCT basis.
    #[must_use]
    pub fn is_prism(&self) -> bool {
        self.cheeks == WindowCheeks::Radial && self.ramp.is_none()
    }
}

/// Declared discretization and budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallWindowOptions {
    /// Cells along the window's axial length.
    pub axial_cells: usize,
    /// Cells across the window's arc.
    pub arc_cells: usize,
    /// Radial layers through the wall.
    pub wall_layers: usize,
    /// Growth ratio of cells away from the window, in `[1, 2]`.
    pub grading: f64,
    /// Downstream port distance from the window centre, in bore radii.
    pub port_radii: f64,
    /// Exterior truncation in bore radii: the outer radius and the axial
    /// half-length.
    pub exterior_radii: f64,
    /// More cells than this along any grid line, or around the bore, are
    /// refused.
    pub max_line_cells: usize,
    /// A coupled face system, or a staircase's dense face map, with more
    /// unknowns than this is refused.
    pub max_unknowns: usize,
    /// Workers for the `theta` modes. Results are bitwise identical at every
    /// count.
    pub workers: usize,
    /// A staircase's conjugate gradients stop at this relative residual.
    pub tolerance: f64,
    /// A staircase needing more conjugate-gradient iterations is refused.
    pub max_iterations: usize,
}

impl Default for WallWindowOptions {
    fn default() -> Self {
        Self {
            axial_cells: 16,
            arc_cells: 24,
            wall_layers: 16,
            grading: 1.1,
            port_radii: 8.0,
            exterior_radii: 10.0,
            max_line_cells: 4096,
            max_unknowns: 4096,
            workers: 1,
            tolerance: 1e-12,
            max_iterations: 100_000,
        }
    }
}

impl WallWindowOptions {
    /// One ladder rung for `window` with `axial_cells` along it, and every
    /// other option at its default.
    ///
    /// The arc cells are about square: within two of
    /// `axial_cells * arc_width / axial_length`, chosen so that the rungs of
    /// [`wall_window_ladder`] (this one and its doubling; for a staircase also
    /// its quadrupling) resolve the requested arc most nearly. The wall layers
    /// are no thicker than an axial cell, and there are at least
    /// `axial_cells` of them.
    #[must_use]
    pub fn rung(window: &WallWindow, axial_cells: usize) -> Self {
        let theta_w = 0.5 * window.arc_width / window.bore_radius;
        let square = (axial_cells as f64 * window.arc_width / window.axial_length)
            .round()
            .max(2.0) as usize;
        let miss = |k: usize| {
            let exact = PI * k as f64 / theta_w;
            (exact.round() - exact).abs() / exact
        };
        let rungs = |k: usize| {
            let doubled = miss(k) + miss(2 * k);
            if window.is_prism() {
                doubled
            } else {
                doubled + miss(4 * k)
            }
        };
        let arc_cells = (square.saturating_sub(2).max(2)..=square + 2)
            .min_by(|&p, &q| {
                rungs(p)
                    .total_cmp(&rungs(q))
                    .then(p.abs_diff(square).cmp(&q.abs_diff(square)))
            })
            .unwrap_or(square);
        let per_cell = window.wall_thickness * axial_cells as f64 / window.axial_length;
        Self {
            axial_cells,
            arc_cells,
            wall_layers: (per_cell.ceil() as usize).max(axial_cells),
            ..Self::default()
        }
    }
}

/// The rungs of one window and their extrapolation.
#[derive(Debug, Clone, PartialEq)]
pub struct WallWindowLadder {
    /// The declared rung first; each later one doubles the cells of the one
    /// before in every direction. Two rungs for a prism, three for a
    /// staircase.
    pub rungs: Vec<WallWindowMass>,
    /// The order the extrapolation used: the re-entrant edges' 4/3 for a
    /// prism, the order observed over the three rungs for a staircase.
    pub order: f64,
    /// `C / Q` extrapolated [1/m].
    pub mass_per_density: f64,
}

/// `2^(4/3) - 1`: Richardson's divisor for one doubling at order 4/3.
const EDGE_RICHARDSON: f64 = 1.519_842_099_789_746_3;

/// A staircase ladder whose observed order lies outside this range is
/// refused. First-order cut cells and the smooth second order of the
/// interior bracket the asymptotic order, which the edges' 4/3 lies between;
/// the margin below admits rungs that approach it from either side, not
/// rungs that have not reached the asymptotic range.
const STAIRCASE_ORDERS: std::ops::RangeInclusive<f64> = 0.5..=2.0;

/// Solve `options` and its doublings, and extrapolate the mass: a prism over
/// two rungs at the re-entrant edges' order 4/3, a staircase over three at
/// the order they show.
///
/// # Errors
/// Those of [`wall_window_mass`] at any rung;
/// [`WallWindowError::Invalid`] when the rungs resolve arcs more than 1e-3
/// apart; [`WallWindowError::Unphysical`] when the masses do not fall toward
/// a limit, or a staircase's observed order lies outside `[0.5, 2]`.
pub fn wall_window_ladder(
    window: &WallWindow,
    options: &WallWindowOptions,
) -> Result<WallWindowLadder> {
    let count = if window.is_prism() { 2 } else { 3 };
    let mut rungs = vec![wall_window_mass(window, options)?];
    let mut cells = *options;
    for _ in 1..count {
        let doubled = |n: usize| {
            n.checked_mul(2).ok_or(WallWindowError::Invalid(
                "cell counts overflow when doubled",
            ))
        };
        cells = WallWindowOptions {
            axial_cells: doubled(cells.axial_cells)?,
            arc_cells: doubled(cells.arc_cells)?,
            wall_layers: doubled(cells.wall_layers)?,
            ..cells
        };
        rungs.push(wall_window_mass(window, &cells)?);
    }
    if rungs
        .iter()
        .any(|r| (r.resolved_arc_width / rungs[0].resolved_arc_width - 1.0).abs() > 1e-3)
    {
        return Err(WallWindowError::Invalid(
            "the rungs resolve different arcs; choose arc cells whose theta cells round alike",
        ));
    }
    let drops: Vec<f64> = rungs
        .windows(2)
        .map(|pair| pair[0].mass_per_density - pair[1].mass_per_density)
        .collect();
    if drops.iter().any(|&d| !(d > 0.0)) {
        return Err(WallWindowError::Unphysical(
            "the ladder does not fall toward its limit; refine the coarse rung",
        ));
    }
    let finest = rungs[count - 1].mass_per_density;
    let (order, mass_per_density) = if count == 2 {
        (4.0 / 3.0, finest - drops[0] / EDGE_RICHARDSON)
    } else {
        let ratio = drops[0] / drops[1];
        let order = det::ln(ratio) / core::f64::consts::LN_2;
        if !STAIRCASE_ORDERS.contains(&order) {
            return Err(WallWindowError::Unphysical(
                "the staircase ladder's observed order lies outside [0.5, 2]; refine the coarse rung",
            ));
        }
        (order, finest - drops[1] / (ratio - 1.0))
    };
    Ok(WallWindowLadder {
        rungs,
        order,
        mass_per_density,
    })
}

/// A solved window mass and its discretization.
#[derive(Debug, Clone, PartialEq)]
pub struct WallWindowMass {
    /// `C / Q` [1/m]: the window branch's acoustic mass per unit density.
    pub mass_per_density: f64,
    /// The same mass as a length of bore, `C S / Q` [m].
    pub bore_length: f64,
    /// The same mass as a length at the resolved window area, `C S_w / Q` [m].
    pub window_length: f64,
    /// Resolved window area on the bore surface, `S_w` [m^2].
    pub window_area: f64,
    /// Open area of the window's outer face [m^2]: each cell open over the
    /// smaller of its fractions at the last wall layer's mid-radius and at
    /// the outer surface (a release face: at the mid-radius), on the outer
    /// surface.
    pub outer_area: f64,
    /// The requested arc as resolved by whole `theta` cells [m].
    pub resolved_arc_width: f64,
    /// Length of the window's outer face along the axis [m]: the window
    /// length plus the ramp's run at the last wall layer's mid-radius.
    pub outer_axial_length: f64,
    /// `Q`: flux into the window's bore face per unit potential [m].
    pub flux: f64,
    /// Flux out of the window's outer face per unit potential [m]. It equals
    /// `flux` to rounding.
    pub outer_flux: f64,
    /// The outer face's open cells on the outer surface, both sides of the
    /// mid-plane, and the flux through each: how `outer_flux` spreads over
    /// the opening. Their fluxes sum to `outer_flux` and their areas to
    /// `outer_area`.
    pub outer_cells: Vec<OuterFaceCell>,
    /// `theta` cells around the bore, `N`.
    pub theta_cells: usize,
    /// Bore grid cells, (radial, axial).
    pub bore_cells: (usize, usize),
    /// Exterior grid cells, (radial, axial); zero with a release face.
    pub exterior_cells: (usize, usize),
    /// Unknowns of the final solve: a prism's coupled face system, a
    /// staircase's folded channel cells and face potentials.
    pub unknowns: usize,
    /// Relative residual of the final solve, recomputed after it.
    pub residual: f64,
    /// Conjugate-gradient iterations; zero for a prism.
    pub iterations: usize,
}

/// The open part of one cell of the window's outer face, on the outer
/// surface `r = a + t`, and the flux through it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OuterFaceCell {
    /// Angle of the open part's centre from the window's mid-plane [rad].
    pub theta: f64,
    /// Angle the open part spans [rad].
    pub arc: f64,
    /// Axial position of the open part's centre from the window centre,
    /// downstream positive [m].
    pub axial: f64,
    /// Axial length of the open part [m].
    pub length: f64,
    /// Flux out through the open part per unit potential [m].
    pub flux: f64,
}

/// Why a window mass was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum WallWindowError {
    /// A declared input is outside its domain.
    Invalid(&'static str),
    /// A declared budget is smaller than the problem needs.
    Budget {
        /// The budgeted quantity.
        what: &'static str,
        /// What the problem needs.
        required: usize,
        /// What was declared.
        maximum: usize,
    },
    /// A factorization failed.
    Factor {
        /// Which factorization.
        stage: &'static str,
        /// The structural failure.
        source: FactorError,
    },
    /// The solve produced a non-positive flux or mass.
    Unphysical(&'static str),
}

impl fmt::Display for WallWindowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "invalid wall window input: {what}"),
            Self::Budget {
                what,
                required,
                maximum,
            } => write!(f, "{what}: {required} exceeds the declared {maximum}"),
            Self::Factor { stage, source } => write!(f, "{stage}: {source}"),
            Self::Unphysical(what) => write!(f, "wall window solve: {what}"),
        }
    }
}

impl std::error::Error for WallWindowError {}

type Result<T> = std::result::Result<T, WallWindowError>;

/// Solve the acoustic mass of `window` under `options`.
///
/// # Errors
/// [`WallWindowError::Invalid`] for a non-physical geometry or option;
/// [`WallWindowError::Budget`] when a declared budget is exceeded;
/// [`WallWindowError::Factor`] when a factorization fails;
/// [`WallWindowError::Unphysical`] when the solved flux or mass is not
/// positive.
#[allow(clippy::too_many_lines)] // one admission, assembly and solve narrative
pub fn wall_window_mass(
    window: &WallWindow,
    options: &WallWindowOptions,
) -> Result<WallWindowMass> {
    let (w, o) = (window, options);
    let positive = |v: f64| v.is_finite() && v > 0.0;
    if !(positive(w.bore_radius)
        && positive(w.wall_thickness)
        && positive(w.arc_width)
        && positive(w.axial_length))
    {
        return Err(WallWindowError::Invalid(
            "bore radius, wall thickness, arc width and axial length must be finite and positive",
        ));
    }
    let (a, t, len) = (w.bore_radius, w.wall_thickness, w.axial_length);
    let half = 0.5 * len;
    let (upstream, through) = match w.upstream {
        WindowUpstream::ClosedEnd { distance } => (distance, false),
        WindowUpstream::Port { distance } => (distance, true),
    };
    if !(upstream.is_finite() && upstream >= half) || (through && !(upstream > half)) {
        return Err(WallWindowError::Invalid(
            "the upstream end must be finite and at least half the window length from its centre, a port beyond it",
        ));
    }
    if o.axial_cells == 0 || o.arc_cells == 0 || o.wall_layers == 0 || o.workers == 0 {
        return Err(WallWindowError::Invalid(
            "cell counts, wall layers and workers must be positive",
        ));
    }
    if !(o.grading >= 1.0 && o.grading <= 2.0) {
        return Err(WallWindowError::Invalid("grading must lie in [1, 2]"));
    }
    let port = o.port_radii * a;
    if !(port.is_finite() && port > half) {
        return Err(WallWindowError::Invalid(
            "the downstream port must lie beyond the window",
        ));
    }
    let free = w.outside == WindowOutside::Free;
    let reach = o.exterior_radii * a;
    if free && !(reach.is_finite() && reach > (a + t).max(half)) {
        return Err(WallWindowError::Invalid(
            "the exterior must reach beyond the wall and the window",
        ));
    }
    let theta_w = 0.5 * w.arc_width / a;
    if !(theta_w < PI) {
        return Err(WallWindowError::Invalid(
            "the window arc must be shorter than the bore's circumference",
        ));
    }
    if w.cheeks == WindowCheeks::Parallel && !(theta_w < 0.5 * PI) {
        return Err(WallWindowError::Invalid(
            "parallel cheeks must meet the bore on less than half its circumference",
        ));
    }
    if w.ramp.is_some_and(|ramp| {
        !(ramp.edge_height >= 0.0
            && ramp.edge_height <= t
            && ramp.angle > 0.0
            && ramp.angle < 0.5 * PI)
    }) {
        return Err(WallWindowError::Invalid(
            "a ramp's edge height must lie in [0, t] and its angle in (0, pi/2)",
        ));
    }
    if !(o.tolerance > 0.0 && o.tolerance < 1.0) || o.max_iterations == 0 {
        return Err(WallWindowError::Invalid(
            "the conjugate-gradient tolerance must lie in (0, 1) and its iteration budget be positive",
        ));
    }
    let around = (PI * o.arc_cells as f64 / theta_w).round();
    if around > o.max_line_cells as f64 {
        return Err(WallWindowError::Budget {
            what: "theta cells around the bore",
            required: around.min(usize::MAX as f64) as usize,
            maximum: o.max_line_cells,
        });
    }
    let n_theta = around as usize;
    if n_theta <= o.arc_cells {
        return Err(WallWindowError::Invalid(
            "the window arc leaves no wall around the bore at this resolution",
        ));
    }
    let (nwx, nwth) = (o.axial_cells, o.arc_cells);
    let dth = 2.0 * PI / n_theta as f64;
    let dxw = len / nwx as f64;
    let layer = t / o.wall_layers as f64;
    let edge_r: Vec<f64> = (0..=o.wall_layers)
        .map(|l| a + t * l as f64 / o.wall_layers as f64)
        .collect();
    // A staircase's cut per wall layer; a prism's layers are the window.
    let layers = if w.is_prism() {
        Vec::new()
    } else {
        staircase_layers(w, &edge_r, nwth, nwx, dth, dxw, o.max_line_cells)?
    };
    // A staircase's outer face against free space: the last layer's cells,
    // each open over the smaller of its own fraction and the cut's at the
    // outer surface, so that the exterior sees the opening the wall has.
    let outer: Vec<f64> = match layers.last() {
        Some(top) if free && w.cheeks == WindowCheeks::Parallel => {
            let half_chord = a * det::sin(0.5 * nwth as f64 * dth);
            let surface = cheek_fractions(det::asin(half_chord / (a + t)), dth, nwth);
            top.arc
                .iter()
                .zip(&surface)
                .map(|(&f, &s)| f.min(s))
                .collect()
        }
        Some(top) => top.arc.clone(),
        None => Vec::new(),
    };
    if !layers.is_empty() && (outer.is_empty() || layers.iter().any(|cut| cut.arc.is_empty())) {
        return Err(WallWindowError::Invalid(
            "the cut closes inside the wall at this resolution",
        ));
    }
    // The outer face's (arc, axial) cells, and its open arc cells and length.
    let (outer_open, outer_arc, outer_axial_length) = match layers.last() {
        Some(top) => (
            (2 * outer.len() - nwth % 2, top.axial.len()),
            (0..outer.len())
                .map(|i| multiplicity(nwth, i) * outer[i])
                .sum::<f64>(),
            dxw * top.axial.iter().sum::<f64>(),
        ),
        None => ((nwth, nwx), nwth as f64, dxw * nwx as f64),
    };
    if free && !(reach > outer_axial_length - half) {
        return Err(WallWindowError::Invalid(
            "the exterior must reach beyond the ramp",
        ));
    }
    let faces = [(nwth, nwx), outer_open].map(|(arc, axial)| FaceBasis::new(arc, axial));
    if w.is_prism() {
        let face = faces[0].size();
        let unknowns = if free { 2 * face } else { face };
        if unknowns > o.max_unknowns {
            return Err(WallWindowError::Budget {
                what: "coupled face unknowns",
                required: unknowns,
                maximum: o.max_unknowns,
            });
        }
    } else {
        let dense = if free {
            faces[0].size().max(faces[1].size())
        } else {
            faces[0].size()
        };
        if dense > o.max_unknowns {
            return Err(WallWindowError::Budget {
                what: "dense face unknowns",
                required: dense,
                maximum: o.max_unknowns,
            });
        }
    }

    let h = dxw.min(a * dth);
    let (ratio, cap) = (o.grading, o.max_line_cells);
    // Uniform cells over the face from the window's upstream edge, the last
    // one ending where the face's opening does, graded beyond; the window's
    // own downstream edge is exactly `half`.
    let axial = |far_up: f64,
                 far_down: f64,
                 pocket: bool,
                 cells: usize,
                 last: f64|
     -> Result<(Vec<f64>, usize)> {
        let mut edges = Vec::new();
        if pocket {
            let mut up = graded(-half, -far_up, dxw, ratio, cap)?;
            up.reverse();
            edges.extend_from_slice(&up[..up.len() - 1]);
        }
        let first = edges.len();
        edges.extend((0..cells).map(|i| -half + len * i as f64 / nwx as f64));
        edges.push(if last < 1.0 {
            edges[edges.len() - 1] + last * dxw
        } else {
            -half + len * cells as f64 / nwx as f64
        });
        let start = if cells == nwx && last == 1.0 {
            half
        } else {
            edges[edges.len() - 1]
        };
        edges.extend_from_slice(&graded(start, far_down, dxw, ratio, cap)?[1..]);
        Ok((edges, first))
    };
    let (bore_x, bore_first) = axial(upstream, port, upstream - half > 1e-6 * dxw, nwx, 1.0)?;
    let mut bore_r = graded(a, 0.0, h, ratio, cap)?;
    bore_r.reverse();
    let bore = Region::new(
        &bore_r,
        &bore_x,
        bore_r.len() - 2,
        bore_first..bore_first + nwx,
        Far::Ports { upstream: through },
    );
    let exterior = if free {
        let last = layers
            .last()
            .map_or(1.0, |top| top.axial[top.axial.len() - 1]);
        let (ext_x, ext_first) = axial(reach, reach, true, outer_open.1, last)?;
        let ext_r = graded(a + t, reach, h, ratio, cap)?;
        Some(Region::new(
            &ext_r,
            &ext_x,
            0,
            ext_first..ext_first + outer_open.1,
            Far::Robin,
        ))
    } else {
        None
    };

    // A staircase's exterior has its own `theta` modes (cut_exterior_map).
    let context = Modes {
        n_theta,
        dth,
        faces: &faces,
        bore: &bore,
        exterior: exterior.as_ref().filter(|_| w.is_prism()),
    };
    let terms = in_parallel(n_theta / 2 + 1, o.workers, |m| context.terms(m))?;

    // Face maps in each face's even DCT basis: potential per unit flux out
    // of the region through each mode, plus the face cells' half-cell drop.
    let face_map = |basis: &FaceBasis,
                    pick: fn(&ModeTerms) -> (&[f64], &[f64]),
                    half_cell: f64,
                    radius: f64| {
        let (even, nx) = (basis.even, basis.axial_cells);
        let face = basis.size();
        let mut map = vec![0.0; face * face];
        for (m, term) in terms.iter().enumerate() {
            let weight = if m == 0 || 2 * m == n_theta { 1.0 } else { 2.0 };
            let g = weight / (n_theta as f64 * dth);
            let (arc, y) = pick(term);
            for pe in 0..even {
                for pf in 0..even {
                    let coefficient = g * arc[pe] * arc[pf];
                    for q in 0..nx {
                        let row = &mut map[(pe * nx + q) * face + pf * nx..][..nx];
                        for (entry, &yv) in row.iter_mut().zip(&y[q * nx..(q + 1) * nx]) {
                            *entry += coefficient * yv;
                        }
                    }
                }
            }
        }
        let diagonal = half_cell / (radius * dth * dxw);
        for i in 0..face {
            map[i * face + i] += diagonal;
        }
        map
    };
    let bore_map = face_map(&faces[0], bore_terms, bore.half_cell, a);
    let exterior_map = match (&exterior, layers.last()) {
        (Some(region), None) => Some(face_map(&faces[1], exterior_terms, region.half_cell, a + t)),
        (Some(region), Some(top)) => Some(cut_exterior_map(
            region,
            &OuterFace {
                basis: &faces[1],
                arc: &outer,
                axial: &top.axial,
                radius: a + t,
            },
            n_theta,
            dth,
            dxw,
            o.workers,
        )?),
        (None, _) => None,
    };

    let solved = if w.is_prism() {
        let (even, face) = (faces[0].even, faces[0].size());
        // Channel: per even DCT mode, flux into the bottom face and out of the
        // top face per unit bottom and top potential, scaled to cell fluxes.
        let mut into_bottom = vec![0.0; face];
        let mut cross = vec![0.0; face];
        let mut out_top = vec![0.0; face];
        let (gb, gt) = (a / (0.5 * layer), (a + t) / (0.5 * layer));
        let scale = dth * dxw;
        for pe in 0..even {
            let mu_p = (2.0 - 2.0 * det::cos(PI * (2 * pe) as f64 / nwth as f64)) / (dth * dth);
            for q in 0..nwx {
                let mu_q = (2.0 - 2.0 * det::cos(PI * q as f64 / nwx as f64)) / (dxw * dxw);
                let (k00, k0l, kll) = channel_inverse_corners(&edge_r, layer, mu_p, mu_q, gb, gt);
                let i = pe * nwx + q;
                into_bottom[i] = scale * (gb - gb * gb * k00);
                cross[i] = -scale * gb * gt * k0l;
                out_top[i] = scale * (gt * gt * kll - gt);
            }
        }

        // [I + N_b B, N_b X; N_e X, I - N_e T] [phi_b; phi_t] = [sqrt(n) e_0; 0],
        // with B, X, T the channel's diagonal maps and N_b, N_e the face maps.
        let size = if free { 2 * face } else { face };
        let mut system = vec![0.0; size * size];
        for r in 0..face {
            for c in 0..face {
                let nb = bore_map[r * face + c];
                system[r * size + c] = nb * into_bottom[c];
                if let Some(ne) = &exterior_map {
                    system[r * size + face + c] = nb * cross[c];
                    system[(face + r) * size + c] = ne[r * face + c] * cross[c];
                    system[(face + r) * size + face + c] = -ne[r * face + c] * out_top[c];
                }
            }
        }
        for i in 0..size {
            system[i * size + i] += 1.0;
        }
        let mut rhs = vec![0.0; size];
        rhs[0] = det::sqrt((nwth * nwx) as f64);
        let factor = lu(&system, size).map_err(|source| WallWindowError::Factor {
            stage: "coupled face system",
            source,
        })?;
        let mut phi = rhs.clone();
        factor.solve(&mut phi);
        let mut misfit = 0.0;
        for (r, row) in system.chunks_exact(size).enumerate() {
            let v: f64 = row.iter().zip(&phi).map(|(m, x)| m * x).sum::<f64>() - rhs[r];
            misfit += v * v;
        }
        let top = |i: usize| if free { phi[face + i] } else { 0.0 };
        let outer_modes: Vec<f64> = (0..face)
            .map(|i| -cross[i] * phi[i] + out_top[i] * top(i))
            .collect();
        Solved {
            flux: rhs[0] * (into_bottom[0] * phi[0] + cross[0] * top(0)),
            outer_flux: rhs[0] * outer_modes[0],
            outer_cells: faces[0].half_cells(&outer_modes),
            unknowns: size,
            residual: det::sqrt(misfit) / rhs[0],
            iterations: 0,
        }
    } else {
        staircase_channel(
            &Staircase {
                layers: &layers,
                outer: &outer,
                edge_r: &edge_r,
                layer,
                dth,
                dxw,
                faces: &faces,
                conductance: [a, a + t].map(|r| r * dth * dxw / (0.5 * layer)),
                resistance: [
                    bore.half_cell / (a * dth * dxw),
                    exterior
                        .as_ref()
                        .map_or(0.0, |e| e.half_cell / ((a + t) * dth * dxw)),
                ],
            },
            &bore_map,
            exterior_map.as_deref(),
            o,
        )?
    };

    let Solved {
        flux,
        outer_flux,
        outer_cells: half_fluxes,
        unknowns,
        residual,
        iterations,
    } = solved;
    // Each half cell's open part, from its open fractions: an arc cell is
    // open from its edge nearer the mid-plane (the odd arc's mid-plane cell
    // about its centre), an axial cell from its upstream edge.
    let (outer_arc_cells, outer_axial) = match layers.last() {
        Some(top) => (outer.clone(), top.axial.clone()),
        None => (vec![1.0; nwth.div_ceil(2)], vec![1.0; nwx]),
    };
    let odd = nwth % 2 == 1;
    let mut outer_cells = Vec::with_capacity(2 * half_fluxes.len());
    for (i, &fi) in outer_arc_cells.iter().enumerate() {
        let theta = if odd && i == 0 {
            0.0
        } else {
            (if odd { i as f64 - 0.5 } else { i as f64 } + 0.5 * fi) * dth
        };
        for (j, &fj) in outer_axial.iter().enumerate() {
            let cell = OuterFaceCell {
                theta,
                arc: fi * dth,
                axial: -half + (j as f64 + 0.5 * fj) * dxw,
                length: fj * dxw,
                flux: half_fluxes[i * outer_axial.len() + j],
            };
            outer_cells.push(cell);
            if !(odd && i == 0) {
                outer_cells.push(OuterFaceCell {
                    theta: -theta,
                    ..cell
                });
            }
        }
    }
    let bore_area = PI * a * a;
    let extrapolated = if through {
        1.0 - flux / (bore_area * (1.0 / port + 1.0 / upstream))
    } else {
        1.0 - flux * port / bore_area
    };
    if !(flux > 0.0 && extrapolated > 0.0) {
        return Err(WallWindowError::Unphysical(
            "the flux or the extrapolated junction potential is not positive",
        ));
    }
    let mass_per_density = extrapolated / flux;
    let resolved_arc_width = nwth as f64 * dth * a;
    let window_area = resolved_arc_width * len;
    Ok(WallWindowMass {
        mass_per_density,
        bore_length: mass_per_density * bore_area,
        window_length: mass_per_density * window_area,
        window_area,
        outer_area: outer_arc * dth * (a + t) * outer_axial_length,
        resolved_arc_width,
        outer_axial_length,
        flux,
        outer_flux,
        outer_cells,
        theta_cells: n_theta,
        bore_cells: (bore.nr, bore.nx),
        exterior_cells: exterior.as_ref().map_or((0, 0), |e| (e.nr, e.nx)),
        unknowns,
        residual,
        iterations,
    })
}

/// The channel's result: what the mass is formed from.
struct Solved {
    flux: f64,
    outer_flux: f64,
    /// Flux out of each of the outer face's half cells, outward from the
    /// mid-plane by axial cells.
    outer_cells: Vec<f64>,
    unknowns: usize,
    residual: f64,
    iterations: usize,
}

/// One face's cells and orthonormal DCT-II bases.
struct FaceBasis {
    arc_cells: usize,
    axial_cells: usize,
    /// The `theta`-even modes: `ceil(arc_cells / 2)`.
    even: usize,
    arc_dct: Vec<f64>,
    axial_dct: Vec<f64>,
}

impl FaceBasis {
    fn new(arc_cells: usize, axial_cells: usize) -> Self {
        Self {
            arc_cells,
            axial_cells,
            even: arc_cells.div_ceil(2),
            arc_dct: dct_rows(arc_cells),
            axial_dct: dct_rows(axial_cells),
        }
    }

    /// Unknowns of the face in its even basis.
    fn size(&self) -> usize {
        self.even * self.axial_cells
    }

    /// The cell values, on the half cells outward from the mid-plane by
    /// axial cells, of the face vector whose even coefficients are `modes`.
    fn half_cells(&self, modes: &[f64]) -> Vec<f64> {
        let (n, nx) = (self.arc_cells, self.axial_cells);
        let half = n.div_ceil(2);
        let mut cells = vec![0.0; half * nx];
        for pe in 0..self.even {
            for q in 0..nx {
                let axial = &self.axial_dct[q * nx..(q + 1) * nx];
                for i in 0..half {
                    let z = modes[pe * nx + q] * self.arc_dct[2 * pe * n + n / 2 + i];
                    for (c, &p) in cells[i * nx..(i + 1) * nx].iter_mut().zip(axial) {
                        *c += z * p;
                    }
                }
            }
        }
        cells
    }

    /// The even arc DCT coefficients of `cos(m theta)` sampled at the arc
    /// cells' centres, with `n` cells around the bore.
    fn arc_coefficients(&self, m: usize, n: f64) -> Vec<f64> {
        let nwth = self.arc_cells;
        let cos_m: Vec<f64> = (0..nwth)
            .map(|k| {
                let centre = k as f64 - 0.5 * (nwth - 1) as f64;
                det::cos(2.0 * PI * m as f64 * centre / n)
            })
            .collect();
        (0..self.even)
            .map(|pe| {
                let row = &self.arc_dct[2 * pe * nwth..(2 * pe + 1) * nwth];
                row.iter().zip(&cos_m).map(|(p, c)| p * c).sum::<f64>()
            })
            .collect()
    }
}

/// What every `theta` mode's terms share.
struct Modes<'a> {
    n_theta: usize,
    dth: f64,
    /// The bore face and the outer face.
    faces: &'a [FaceBasis; 2],
    bore: &'a Region,
    exterior: Option<&'a Region>,
}

impl Modes<'_> {
    /// Mode `m`: the discrete `theta` Laplacian's eigenvalue `mu`, each
    /// face's arc coefficients of `cos(m theta)`, and each region's face
    /// response projected onto its face's axial DCT basis.
    fn terms(&self, m: usize) -> Result<ModeTerms> {
        let n = self.n_theta as f64;
        let mu = (2.0 - 2.0 * det::cos(2.0 * PI * m as f64 / n)) / (self.dth * self.dth);
        let [bore_face, outer_face] = self.faces;
        let response = |region: &Region, face: &FaceBasis, stage: &'static str| {
            region
                .face_response(mu)
                .map(|y| congruence(&face.axial_dct, &y, face.axial_cells))
                .map_err(|source| WallWindowError::Factor { stage, source })
        };
        let bore = response(self.bore, bore_face, "bore mode elimination")?;
        let exterior = match self.exterior {
            Some(region) => response(region, outer_face, "exterior mode elimination")?,
            None => Vec::new(),
        };
        Ok(ModeTerms {
            bore_arc: bore_face.arc_coefficients(m, n),
            outer_arc: outer_face.arc_coefficients(m, n),
            bore,
            exterior,
        })
    }
}

/// One `theta` mode's contribution: each face's even arc DCT coefficients of
/// `cos(m theta)`, and each region's face response in its face's axial DCT
/// basis (axial cells squared, row-major).
struct ModeTerms {
    bore_arc: Vec<f64>,
    outer_arc: Vec<f64>,
    bore: Vec<f64>,
    exterior: Vec<f64>,
}

fn bore_terms(term: &ModeTerms) -> (&[f64], &[f64]) {
    (&term.bore_arc, &term.bore)
}

fn exterior_terms(term: &ModeTerms) -> (&[f64], &[f64]) {
    (&term.outer_arc, &term.exterior)
}

/// `f` over `0..count` in order, on up to `workers` threads.
fn in_parallel<T: Send>(
    count: usize,
    workers: usize,
    f: impl Fn(usize) -> Result<T> + Sync,
) -> Result<Vec<T>> {
    let per = count.div_ceil(workers.min(count)).max(1);
    let chunks: Vec<Result<Vec<T>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..count)
            .step_by(per)
            .map(|start| {
                let f = &f;
                scope.spawn(move || (start..(start + per).min(count)).map(f).collect())
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("wall window mode worker panicked"))
            .collect()
    });
    let mut out = Vec::with_capacity(count);
    for chunk in chunks {
        out.extend(chunk?);
    }
    Ok(out)
}

/// A staircase's outer face as the exterior sees it.
struct OuterFace<'a> {
    basis: &'a FaceBasis,
    /// Open fraction of each half arc cell, outward from the mid-plane.
    arc: &'a [f64],
    /// Open fraction of each axial cell.
    axial: &'a [f64],
    /// The outer surface's radius [m].
    radius: f64,
}

/// The exterior's face map under a staircase's outer face, in the face's
/// even basis. The face's last open arc cell can be cut, so the exterior has
/// its own `theta` grid: the bore's cells folded onto `[0, pi]`, with that
/// cell split where its opening ends. Its modes solve `K v = mu M v`, `M` the
/// cell widths and `K` the links of one over the centre spacing; with
/// `v^T M v = 1` over the half circle, `v / sqrt(2)` is the mode normalized
/// over the whole circle. The exterior's axial grid already ends the face's
/// last cell where its opening does. Each face cell's half-cell drop covers
/// its open part only.
///
/// # Errors
/// [`WallWindowError::Invalid`] when the folded grid exceeds the dense
/// eigensolver's workspace; [`WallWindowError::Factor`] when a mode's
/// elimination fails.
fn cut_exterior_map(
    region: &Region,
    face: &OuterFace<'_>,
    n_theta: usize,
    dth: f64,
    dxw: f64,
    workers: usize,
) -> Result<Vec<f64>> {
    let OuterFace {
        basis,
        arc,
        axial,
        radius,
    } = *face;
    let (even, nx, narc) = (basis.even, basis.axial_cells, basis.arc_cells);
    let edge = arc.len() - 1;
    let (widths, centres) = folded_cells(n_theta, dth, narc % 2 == 1, edge, arc[edge]);
    let n = widths.len();
    if fs_la::eigen::admit_jacobi_eigh(n).is_err() {
        return Err(WallWindowError::Invalid(
            "the exterior's theta cells exceed the dense eigensolver's workspace",
        ));
    }
    // M^(-1/2) K M^(-1/2), whose eigenvectors are M^(1/2) v.
    let scale: Vec<f64> = widths.iter().map(|&w| 1.0 / det::sqrt(w)).collect();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n - 1 {
        let g = 1.0 / (centres[i + 1] - centres[i]);
        matrix[i * n + i] += g * scale[i] * scale[i];
        matrix[(i + 1) * n + i + 1] += g * scale[i + 1] * scale[i + 1];
        let c = -g * scale[i] * scale[i + 1];
        matrix[i * n + i + 1] = c;
        matrix[(i + 1) * n + i] = c;
    }
    let (values, vectors) = fs_la::eigen::householder_eigh(&matrix, n);
    // Each face arc cell's folded cell, outward from the mid-plane.
    let folded: Vec<usize> = (0..narc).map(|k| (2 * k).abs_diff(narc - 1) / 2).collect();
    let coefficients: Vec<Vec<f64>> = (0..n)
        .map(|mode| {
            (0..even)
                .map(|pe| {
                    let row = &basis.arc_dct[2 * pe * narc..(2 * pe + 1) * narc];
                    let v: f64 = row
                        .iter()
                        .zip(&folded)
                        .map(|(p, &i)| p * scale[i] * vectors[i * n + mode])
                        .sum();
                    v * core::f64::consts::FRAC_1_SQRT_2
                })
                .collect()
        })
        .collect();
    let responses = in_parallel(n, workers, |mode| {
        region
            .face_response(values[mode].max(0.0))
            .map(|y| congruence(&basis.axial_dct, &y, nx))
            .map_err(|source| WallWindowError::Factor {
                stage: "exterior mode elimination",
                source,
            })
    })?;
    let size = basis.size();
    let mut map = vec![0.0; size * size];
    for (c, y) in coefficients.iter().zip(&responses) {
        for pe in 0..even {
            for pf in 0..even {
                let coefficient = c[pe] * c[pf];
                for q in 0..nx {
                    let row = &mut map[(pe * nx + q) * size + pf * nx..][..nx];
                    for (entry, &yv) in row.iter_mut().zip(&y[q * nx..(q + 1) * nx]) {
                        *entry += coefficient * yv;
                    }
                }
            }
        }
    }
    // The half-cell drops, over each cell's open width.
    let arc_widths: Vec<f64> = folded.iter().map(|&i| arc[i] * dth).collect();
    let axial_widths: Vec<f64> = axial.iter().map(|&f| f * dxw).collect();
    let arc_drop = inverse_widths(&basis.arc_dct, 2, even, &arc_widths);
    let axial_drop = inverse_widths(&basis.axial_dct, 1, nx, &axial_widths);
    let drop = region.half_cell / radius;
    for pe in 0..even {
        for pf in 0..even {
            let arc_term = drop * arc_drop[pe * even + pf];
            for q in 0..nx {
                let row = &mut map[(pe * nx + q) * size + pf * nx..][..nx];
                for (entry, &x) in row.iter_mut().zip(&axial_drop[q * nx..(q + 1) * nx]) {
                    *entry += arc_term * x;
                }
            }
        }
    }
    Ok(map)
}

/// `P diag(1 / widths) P^T` for the `count` DCT rows `0, stride, 2 stride,
/// ...` of `rows`, which has one column per width. Row-major `count x count`.
fn inverse_widths(rows: &[f64], stride: usize, count: usize, widths: &[f64]) -> Vec<f64> {
    let cells = widths.len();
    let mut out = vec![0.0; count * count];
    for p in 0..count {
        let rp = &rows[p * stride * cells..(p * stride + 1) * cells];
        for q in 0..count {
            let rq = &rows[q * stride * cells..(q * stride + 1) * cells];
            out[p * count + q] = rp
                .iter()
                .zip(rq)
                .zip(widths)
                .map(|((x, y), w)| x * y / w)
                .sum();
        }
    }
    out
}

/// The exterior's `theta` cells folded onto `[0, pi]`, as widths and centres
/// [rad]. They are the bore's cells, centred on the mid-plane for an odd arc
/// and beside it for an even one; a cell centred on `0` or `pi` is halved and
/// keeps its centre. The face's edge cell `edge` is split where its open
/// fraction `open` ends.
fn folded_cells(
    n_theta: usize,
    dth: f64,
    odd: bool,
    edge: usize,
    open: f64,
) -> (Vec<f64>, Vec<f64>) {
    // Cell j is centred at twice/2 dth, twice = 2j for an odd arc, 2j + 1 for
    // an even one.
    let mut cells: Vec<(f64, f64)> = (0..)
        .map(|j: usize| 2 * j + usize::from(!odd))
        .take_while(|&twice| twice <= n_theta)
        .map(|twice| {
            let width = if twice == 0 || twice == n_theta {
                0.5 * dth
            } else {
                dth
            };
            (width, 0.5 * twice as f64 * dth)
        })
        .collect();
    if open < 1.0 {
        let (width, centre) = cells[edge];
        let parts = if odd && edge == 0 {
            // The mid-plane cell, open over |theta| <= open dth / 2.
            [
                (0.5 * open * dth, 0.0),
                (0.5 * (1.0 - open) * dth, 0.25 * (1.0 + open) * dth),
            ]
        } else {
            let lower = centre - 0.5 * width;
            [
                (open * dth, lower + 0.5 * open * dth),
                (
                    (1.0 - open) * dth,
                    lower + (open + 0.5 * (1.0 - open)) * dth,
                ),
            ]
        };
        cells.splice(edge..=edge, parts);
    }
    cells.into_iter().unzip()
}

/// One wall layer of a staircase where the cut crosses its mid-radius: the
/// open fraction of each half arc cell, outward from the mid-plane, and of
/// each axial cell, downstream from the window's upstream edge. Cells the cut
/// misses are left out.
struct Cut {
    arc: Vec<f64>,
    axial: Vec<f64>,
}

/// A staircase channel's cut per wall layer.
///
/// # Errors
/// [`WallWindowError::Budget`] when a layer's axial cells exceed
/// `max_line_cells`.
fn staircase_layers(
    window: &WallWindow,
    edge_r: &[f64],
    arc_cells: usize,
    axial_cells: usize,
    dth: f64,
    dxw: f64,
    max_line_cells: usize,
) -> Result<Vec<Cut>> {
    let a = window.bore_radius;
    let half_chord = a * det::sin(0.5 * arc_cells as f64 * dth);
    let ramp = window.ramp.map(|ramp| {
        let (sin, cos) = det::sin_cos(ramp.angle);
        (ramp.edge_height, cos / sin)
    });
    let mut layers = Vec::with_capacity(edge_r.len() - 1);
    for pair in edge_r.windows(2) {
        let mid = 0.5 * (pair[0] + pair[1]);
        let arc = match window.cheeks {
            WindowCheeks::Radial => vec![1.0; arc_cells.div_ceil(2)],
            WindowCheeks::Parallel => cheek_fractions(det::asin(half_chord / mid), dth, arc_cells),
        };
        // The ramp's run past the window's downstream edge, in axial cells.
        let run = ramp.map_or(0.0, |(edge, cot)| (mid - a - edge).max(0.0) * cot / dxw);
        let run = if (run - run.round()).abs() < CUT_SNAP {
            run.round()
        } else {
            run
        };
        let cells = axial_cells as f64 + run.ceil();
        if !(cells <= max_line_cells as f64) {
            return Err(WallWindowError::Budget {
                what: "axial cells over the outer face",
                required: cells.min(usize::MAX as f64) as usize,
                maximum: max_line_cells,
            });
        }
        let mut axial = vec![1.0; axial_cells];
        axial.extend((0..run.ceil() as usize).map(|j| (run - j as f64).min(1.0)));
        layers.push(Cut { arc, axial });
    }
    Ok(layers)
}

/// A cut within this fraction of a cell edge lies on the edge. The sliver
/// cell it would leave carries nothing, but its near-zero width would
/// ill-condition the operators, the exterior's `theta` eigenproblem above all.
const CUT_SNAP: f64 = 1e-9;

/// `f` moved onto a cell edge when within [`CUT_SNAP`] of one.
fn snap(f: f64) -> f64 {
    if f > 1.0 - CUT_SNAP {
        1.0
    } else if f < CUT_SNAP {
        0.0
    } else {
        f
    }
}

/// Open fraction of each half arc cell, outward from the mid-plane, within
/// the angle `phi` that parallel cheeks subtend at some radius; the cells they
/// close are left out.
fn cheek_fractions(phi: f64, dth: f64, arc_cells: usize) -> Vec<f64> {
    let odd = arc_cells % 2 == 1;
    (0..arc_cells.div_ceil(2))
        .map(|i| {
            snap(if odd && i == 0 {
                (2.0 * phi / dth).min(1.0)
            } else {
                let lower = if odd { i as f64 - 0.5 } else { i as f64 } * dth;
                ((phi - lower) / dth).clamp(0.0, 1.0)
            })
        })
        .take_while(|&f| f > 0.0)
        .collect()
}

/// A staircase channel's geometry and grid.
struct Staircase<'a> {
    /// The cut of each layer, bore to outer surface.
    layers: &'a [Cut],
    /// Open fraction of the outer face's half arc cells: against free space
    /// the smaller of the last layer's and the cut's at the outer surface, a
    /// release face the last layer's.
    outer: &'a [f64],
    /// Layer edge radii, bore to outer surface.
    edge_r: &'a [f64],
    /// Layer thickness [m].
    layer: f64,
    dth: f64,
    dxw: f64,
    /// The bore face (the window's cells) and the outer face (the last
    /// layer's).
    faces: &'a [FaceBasis; 2],
    /// Conductance from an open face cell to its channel cell's centre, bore
    /// and outer face [m].
    conductance: [f64; 2],
    /// The half-cell resistance each face map carries on its diagonal, bore
    /// and exterior [1/m].
    resistance: [f64; 2],
}

/// The folded unknowns: each layer's half arc cells by its axial cells, then
/// the bore face's potentials and, with a free exterior, the outer face's. A
/// half cell stands for itself and its mirror image, except the mid-plane
/// cell of an odd arc.
struct Channel {
    /// First potential of the bore face and of the outer face.
    faces: [usize; 2],
    /// Folded conductances `(unknown, unknown, g)` [m].
    links: Vec<(usize, usize, f64)>,
    /// Conductances from top-layer cells to a pressure-release outer face,
    /// `(cell, g)` [m].
    ground: Vec<(usize, f64)>,
    /// Each unknown's summed conductances [m].
    diagonal: Vec<f64>,
}

/// Mirror images a half cell stands for.
fn multiplicity(arc_cells: usize, half_index: usize) -> f64 {
    if arc_cells % 2 == 1 && half_index == 0 {
        1.0
    } else {
        2.0
    }
}

impl Channel {
    fn new(s: &Staircase<'_>, free: bool) -> Self {
        let (layers, layer) = (s.layers.len(), s.layer);
        let mut start = vec![0];
        for cut in s.layers {
            start.push(start[start.len() - 1] + cut.arc.len() * cut.axial.len());
        }
        let [bore_face, outer_face] = s.faces;
        let faces = [start[layers], start[layers] + bore_face.size()];
        let size = if free {
            faces[1] + outer_face.size()
        } else {
            faces[1]
        };
        let m = |i: usize| multiplicity(bore_face.arc_cells, i);
        let mut links = Vec::new();
        for (l, cut) in s.layers.iter().enumerate() {
            let mid = 0.5 * (s.edge_r[l] + s.edge_r[l + 1]);
            let g_theta = layer * s.dxw / (mid * s.dth);
            let g_axial = mid * layer * s.dth / s.dxw;
            let g_radial = s.edge_r[l + 1] * s.dth * s.dxw / layer;
            let nx = cut.axial.len();
            for (i, &fi) in cut.arc.iter().enumerate() {
                for (j, &fj) in cut.axial.iter().enumerate() {
                    let c = start[l] + i * nx + j;
                    // A link into a cut cell ends at its open part's
                    // centroid. The mid-plane link of an even arc joins a
                    // cell to its own image and carries nothing; the rest
                    // pair up.
                    if let Some(&next) = cut.arc.get(i + 1) {
                        links.push((c, c + nx, 2.0 * fj * g_theta / (0.5 * (1.0 + next))));
                    }
                    if let Some(&next) = cut.axial.get(j + 1) {
                        links.push((c, c + 1, m(i) * fi * g_axial / (0.5 * (1.0 + next))));
                    }
                    if let Some(up) = s.layers.get(l + 1)
                        && let (Some(&ui), Some(&uj)) = (up.arc.get(i), up.axial.get(j))
                    {
                        let open = fi.min(ui) * fj.min(uj);
                        let above = start[l + 1] + i * up.axial.len() + j;
                        links.push((c, above, m(i) * open * g_radial));
                    }
                }
            }
        }
        // Each face cell links to its channel cell over that cell's open
        // fraction; a release face holds its potential at zero.
        let bottom = &s.layers[0];
        let nwx = bore_face.axial_cells;
        for (i, &fi) in bottom.arc.iter().enumerate().take(bore_face.even) {
            for (j, &fj) in bottom.axial.iter().enumerate().take(nwx) {
                let c = i * bottom.axial.len() + j;
                links.push((c, faces[0] + i * nwx + j, m(i) * fi * fj * s.conductance[0]));
            }
        }
        let top = &s.layers[layers - 1];
        let mut ground = Vec::new();
        for (i, &fi) in s.outer.iter().enumerate() {
            for (j, &fj) in top.axial.iter().enumerate() {
                let k = i * top.axial.len() + j;
                let g = m(i) * fi * fj * s.conductance[1];
                if free {
                    links.push((start[layers - 1] + k, faces[1] + k, g));
                } else {
                    ground.push((start[layers - 1] + k, g));
                }
            }
        }
        let mut diagonal = vec![0.0; size];
        for &(p, q, g) in &links {
            diagonal[p] += g;
            diagonal[q] += g;
        }
        for &(c, g) in &ground {
            diagonal[c] += g;
        }
        Self {
            faces,
            links,
            ground,
            diagonal,
        }
    }

    /// `y = K u` over the links.
    fn apply(&self, u: &[f64], y: &mut [f64]) {
        for ((yv, &uv), &d) in y.iter_mut().zip(u).zip(&self.diagonal) {
            *yv = d * uv;
        }
        for &(p, q, g) in &self.links {
            y[p] -= g * u[q];
            y[q] -= g * u[p];
        }
    }
}

/// One face's potentials among the unknowns and their even DCT transform.
struct FaceLink<'a> {
    basis: &'a FaceBasis,
    /// First unknown of the face.
    start: usize,
    /// Half arc cells.
    half: usize,
    /// `m_i P[2 pe][k_i]`, `even x half`: the even arc rows on the half
    /// cells, times the cells' multiplicities.
    arc: Vec<f64>,
}

impl<'a> FaceLink<'a> {
    fn new(basis: &'a FaceBasis, start: usize) -> Self {
        let (n, half) = (basis.arc_cells, basis.arc_cells.div_ceil(2));
        let mut arc = vec![0.0; basis.even * half];
        for pe in 0..basis.even {
            for i in 0..half {
                arc[pe * half + i] = multiplicity(n, i) * basis.arc_dct[2 * pe * n + n / 2 + i];
            }
        }
        Self {
            basis,
            start,
            half,
            arc,
        }
    }

    /// The face cells' potentials in the face's even basis.
    fn to_modes(&self, u: &[f64]) -> Vec<f64> {
        let (even, nx) = (self.basis.even, self.basis.axial_cells);
        let cells = &u[self.start..self.start + self.half * nx];
        let mut arc = vec![0.0; even * nx];
        for pe in 0..even {
            let row = &mut arc[pe * nx..(pe + 1) * nx];
            for (i, &f) in self.arc[pe * self.half..(pe + 1) * self.half]
                .iter()
                .enumerate()
            {
                for (o, &v) in row.iter_mut().zip(&cells[i * nx..(i + 1) * nx]) {
                    *o += f * v;
                }
            }
        }
        let mut modes = vec![0.0; even * nx];
        for pe in 0..even {
            for q in 0..nx {
                modes[pe * nx + q] = self.basis.axial_dct[q * nx..(q + 1) * nx]
                    .iter()
                    .zip(&arc[pe * nx..(pe + 1) * nx])
                    .map(|(p, v)| p * v)
                    .sum();
            }
        }
        modes
    }

    /// Add the folded cell values of `modes` (the transpose of
    /// [`Self::to_modes`]) to the face cells of `y`.
    fn add_from_modes(&self, modes: &[f64], y: &mut [f64]) {
        let (even, nx) = (self.basis.even, self.basis.axial_cells);
        let mut arc = vec![0.0; even * nx];
        for pe in 0..even {
            let row = &mut arc[pe * nx..(pe + 1) * nx];
            for q in 0..nx {
                let z = modes[pe * nx + q];
                for (o, &p) in row
                    .iter_mut()
                    .zip(&self.basis.axial_dct[q * nx..(q + 1) * nx])
                {
                    *o += p * z;
                }
            }
        }
        let cells = &mut y[self.start..self.start + self.half * nx];
        for pe in 0..even {
            let row = &arc[pe * nx..(pe + 1) * nx];
            for (i, &f) in self.arc[pe * self.half..(pe + 1) * self.half]
                .iter()
                .enumerate()
            {
                for (o, &v) in cells[i * nx..(i + 1) * nx].iter_mut().zip(row) {
                    *o += f * v;
                }
            }
        }
    }
}

/// A region's face map `N` in its face's even basis, as the lower Cholesky
/// factor of its symmetric part (row-major). Applied inverted, it is the
/// region's Dirichlet-to-Neumann map on the face.
struct FaceLoad {
    factor: Vec<f64>,
    size: usize,
}

impl FaceLoad {
    fn new(map: &[f64], size: usize, stage: &'static str) -> Result<Self> {
        let mut matrix = vec![0.0; size * size];
        for r in 0..size {
            for c in 0..=r {
                let v = 0.5 * (map[r * size + c] + map[c * size + r]);
                matrix[r * size + c] = v;
                matrix[c * size + r] = v;
            }
        }
        let factor =
            cholesky(&matrix, size).map_err(|source| WallWindowError::Factor { stage, source })?;
        for r in 0..size {
            for c in 0..=r {
                matrix[r * size + c] = factor.l(r, c);
            }
        }
        Ok(Self {
            factor: matrix,
            size,
        })
    }

    /// `N^-1 v`.
    fn apply(&self, v: &[f64]) -> Vec<f64> {
        let (n, factor) = (self.size, &self.factor);
        let mut x = v.to_vec();
        for i in 0..n {
            let row = &factor[i * n..i * n + i];
            let s: f64 = row.iter().zip(&x[..i]).map(|(l, y)| l * y).sum();
            x[i] = (x[i] - s) / factor[i * n + i];
        }
        // Back substitution by rows of L, so each pass reads one row.
        for i in (0..n).rev() {
            x[i] /= factor[i * n + i];
            let xi = x[i];
            for (y, &l) in x[..i].iter_mut().zip(&factor[i * n..i * n + i]) {
                *y -= l * xi;
            }
        }
        x
    }
}

/// Solve a staircase channel with its face potentials by Jacobi-
/// preconditioned conjugate gradients over the folded unknowns.
fn staircase_channel(
    s: &Staircase<'_>,
    bore_map: &[f64],
    exterior_map: Option<&[f64]>,
    options: &WallWindowOptions,
) -> Result<Solved> {
    let channel = Channel::new(s, exterior_map.is_some());
    let n = channel.diagonal.len();
    let [bore_face, outer_face] = s.faces;
    // Each face with its half-cell resistance and, if its cells are cut, the
    // open fractions of its half arc and axial cells.
    let mut faces = vec![(
        FaceLink::new(bore_face, channel.faces[0]),
        FaceLoad::new(bore_map, bore_face.size(), "bore face map")?,
        s.resistance[0],
        None,
    )];
    if let Some(map) = exterior_map {
        let top = &s.layers[s.layers.len() - 1];
        faces.push((
            FaceLink::new(outer_face, channel.faces[1]),
            FaceLoad::new(map, outer_face.size(), "exterior face map")?,
            s.resistance[1],
            Some((s.outer, top.axial.as_slice())),
        ));
    }
    let apply = |u: &[f64], y: &mut [f64]| {
        channel.apply(u, y);
        for (link, load, ..) in &faces {
            link.add_from_modes(&load.apply(&link.to_modes(u)), y);
        }
    };
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p * q).sum::<f64>();

    // The ports hold the bore at potential 1: sqrt(n) times the face's
    // constant mode.
    let bore_cells = bore_face.arc_cells * bore_face.axial_cells;
    let mut uniform = vec![0.0; bore_face.size()];
    uniform[0] = det::sqrt(bore_cells as f64);
    let mut rhs = vec![0.0; n];
    let (bore, bore_load, ..) = &faces[0];
    bore.add_from_modes(&bore_load.apply(&uniform), &mut rhs);
    // A face map is at least its cells' half-cell drops, each over the cell's
    // open part, so its inverse's diagonal on a half cell is at most the
    // cell's multiplicity over its own drop.
    let mut jacobi = channel.diagonal.clone();
    for (link, _, resistance, open) in &faces {
        let nx = link.basis.axial_cells;
        for i in 0..link.half {
            let g = multiplicity(link.basis.arc_cells, i) / resistance;
            for (j, v) in jacobi[link.start + i * nx..link.start + (i + 1) * nx]
                .iter_mut()
                .enumerate()
            {
                *v += open.map_or(g, |(arc, axial)| g * arc[i] * axial[j]);
            }
        }
    }

    let norm = det::sqrt(dot(&rhs, &rhs));
    let mut u = vec![0.0; n];
    let mut r = rhs.clone();
    let mut z: Vec<f64> = r.iter().zip(&jacobi).map(|(x, d)| x / d).collect();
    let mut p = z.clone();
    let mut q = vec![0.0; n];
    let mut rz = dot(&r, &z);
    let mut iterations = 0;
    while det::sqrt(dot(&r, &r)) > options.tolerance * norm {
        if iterations == options.max_iterations {
            return Err(WallWindowError::Budget {
                what: "conjugate-gradient iterations",
                required: iterations + 1,
                maximum: options.max_iterations,
            });
        }
        apply(&p, &mut q);
        let step = rz / dot(&p, &q);
        for ((uv, rv), (&pv, &qv)) in u.iter_mut().zip(&mut r).zip(p.iter().zip(&q)) {
            *uv += step * pv;
            *rv -= step * qv;
        }
        for ((zv, &rv), &d) in z.iter_mut().zip(&r).zip(&jacobi) {
            *zv = rv / d;
        }
        let next = dot(&r, &z);
        let beta = next / rz;
        for (pv, &zv) in p.iter_mut().zip(&z) {
            *pv = zv + beta * *pv;
        }
        rz = next;
        iterations += 1;
    }
    apply(&u, &mut q);
    let misfit: f64 = q.iter().zip(&rhs).map(|(y, b)| (y - b) * (y - b)).sum();

    let mut drop = uniform;
    for (d, w) in drop.iter_mut().zip(bore.to_modes(&u)) {
        *d -= w;
    }
    let (outer_flux, outer_cells) = match faces.get(1) {
        Some((top, load, ..)) => {
            let outer_cells = outer_face.arc_cells * outer_face.axial_cells;
            let modes = load.apply(&top.to_modes(&u));
            (
                det::sqrt(outer_cells as f64) * modes[0],
                outer_face.half_cells(&modes),
            )
        }
        None => {
            // A ground link carries its half cell and its mirror image.
            let nx = outer_face.axial_cells;
            let cells = channel.ground.iter().enumerate().map(|(k, &(c, g))| {
                g * u[c] / multiplicity(outer_face.arc_cells, k / nx)
            });
            (
                channel.ground.iter().map(|&(c, g)| g * u[c]).sum(),
                cells.collect(),
            )
        }
    };
    Ok(Solved {
        flux: det::sqrt(bore_cells as f64) * bore_load.apply(&drop)[0],
        outer_flux,
        outer_cells,
        unknowns: n,
        residual: det::sqrt(misfit) / norm,
        iterations,
    })
}

/// A region's far boundaries.
#[derive(Clone, Copy)]
enum Far {
    /// The bore: a port (potential 0 in the response) downstream, and
    /// upstream as well when `upstream`, else a rigid end.
    Ports { upstream: bool },
    /// The exterior: a monopole Robin condition on the outer cylinder and
    /// both axial ends.
    Robin,
}

/// One `theta`-invariant region's finite-volume operator, negated so that
/// every `theta` mode's matrix is positive definite. Axial block `j` is
/// tridiagonal over the radial cells, with diagonal
/// `diag + mu theta_weight` and off-diagonals `-off`. Blocks `j` and `j + 1`
/// couple through `-axial`.
struct Region {
    nr: usize,
    nx: usize,
    diag: Vec<f64>,
    theta_weight: Vec<f64>,
    off: Vec<f64>,
    axial: Vec<f64>,
    /// Radial index of the window face's cells.
    row: usize,
    /// Axial indices of the window face's cells.
    columns: Range<usize>,
    /// Half the radial width of the face cells [m].
    half_cell: f64,
}

impl Region {
    fn new(r_edges: &[f64], x_edges: &[f64], row: usize, columns: Range<usize>, far: Far) -> Self {
        let centre =
            |e: &[f64]| -> Vec<f64> { e.windows(2).map(|p| 0.5 * (p[0] + p[1])).collect() };
        let width = |e: &[f64]| -> Vec<f64> { e.windows(2).map(|p| p[1] - p[0]).collect() };
        let (rc, dr, xc, dx) = (
            centre(r_edges),
            width(r_edges),
            centre(x_edges),
            width(x_edges),
        );
        let (nr, nx) = (rc.len(), xc.len());
        let radial: Vec<f64> = (0..nr - 1)
            .map(|i| r_edges[i + 1] / (rc[i + 1] - rc[i]))
            .collect();
        let mut diag = vec![0.0; nx * nr];
        let mut theta_weight = vec![0.0; nx * nr];
        let mut off = vec![0.0; nx * (nr - 1)];
        let mut axial = vec![0.0; (nx - 1) * nr];
        for j in 0..nx {
            for i in 0..nr {
                theta_weight[j * nr + i] = dr[i] * dx[j] / rc[i];
            }
            for i in 0..nr - 1 {
                let g = radial[i] * dx[j];
                off[j * (nr - 1) + i] = g;
                diag[j * nr + i] += g;
                diag[j * nr + i + 1] += g;
            }
        }
        for j in 0..nx - 1 {
            for i in 0..nr {
                let g = rc[i] * dr[i] / (xc[j + 1] - xc[j]);
                axial[j * nr + i] = g;
                diag[j * nr + i] += g;
                diag[(j + 1) * nr + i] += g;
            }
        }
        match far {
            Far::Ports { upstream } => {
                for i in 0..nr {
                    diag[(nx - 1) * nr + i] += rc[i] * dr[i] / (0.5 * dx[nx - 1]);
                    if upstream {
                        diag[i] += rc[i] * dr[i] / (0.5 * dx[0]);
                    }
                }
            }
            Far::Robin => {
                // d(phi)/dn = -kappa phi for a monopole at the origin, with
                // the boundary potential half a cell beyond the cell centre.
                let outer = r_edges[nr];
                let half_r = 0.5 * dr[nr - 1];
                for j in 0..nx {
                    let kappa = outer / (outer * outer + xc[j] * xc[j]);
                    diag[j * nr + nr - 1] += outer * dx[j] * kappa / (1.0 + kappa * half_r);
                }
                for (j, end) in [(0, -x_edges[0]), (nx - 1, x_edges[nx])] {
                    let half_x = 0.5 * dx[j];
                    for i in 0..nr {
                        let kappa = end / (rc[i] * rc[i] + end * end);
                        diag[j * nr + i] += rc[i] * dr[i] * kappa / (1.0 + kappa * half_x);
                    }
                }
            }
        }
        Self {
            nr,
            nx,
            diag,
            theta_weight,
            off,
            axial,
            row,
            half_cell: 0.5 * dr[row],
            columns,
        }
    }

    /// Potentials on the window face cells for a unit flux out of the region
    /// through each of them, in `theta` mode `mu`. Row-major, `face x face`
    /// over the face's axial cells.
    fn face_response(&self, mu: f64) -> std::result::Result<Vec<f64>, FactorError> {
        let (nr, nx, row) = (self.nr, self.nx, self.row);
        let nw = self.columns.len();
        let c0 = self.columns.start;
        let block = nr * nr;
        let mut inverses = vec![0.0; nx * block];
        let mut z = vec![0.0; nx * nr * nw];
        let mut s = vec![0.0; block];
        let mut unit = vec![0.0; nr];
        for j in 0..nx {
            s.fill(0.0);
            for i in 0..nr {
                s[i * nr + i] = self.diag[j * nr + i] + mu * self.theta_weight[j * nr + i];
            }
            for i in 0..nr - 1 {
                let v = -self.off[j * (nr - 1) + i];
                s[i * nr + i + 1] = v;
                s[(i + 1) * nr + i] = v;
            }
            let (done, rest) = z.split_at_mut(j * nr * nw);
            let zj = &mut rest[..nr * nw];
            if self.columns.contains(&j) {
                zj[row * nw + (j - c0)] = 1.0;
            }
            if j > 0 {
                // Block elimination: S_j = B_j - A S_{j-1}^{-1} A and
                // z_j = e_j + A S_{j-1}^{-1} z_{j-1}, A = diag(axial_{j-1}).
                let prev = &inverses[(j - 1) * block..j * block];
                let ax = &self.axial[(j - 1) * nr..j * nr];
                let zp = &done[(j - 1) * nr * nw..];
                for p in 0..nr {
                    for q in 0..nr {
                        s[p * nr + q] -= ax[p] * prev[p * nr + q] * ax[q];
                    }
                    let out = &mut zj[p * nw..(p + 1) * nw];
                    for q in 0..nr {
                        let f = ax[p] * prev[p * nr + q];
                        if f != 0.0 {
                            for (o, &v) in out.iter_mut().zip(&zp[q * nw..(q + 1) * nw]) {
                                *o += f * v;
                            }
                        }
                    }
                }
            }
            let factor = cholesky(&s, nr)?;
            let inverse = &mut inverses[j * block..(j + 1) * block];
            for c in 0..nr {
                unit.fill(0.0);
                unit[c] = 1.0;
                factor.solve(&mut unit);
                for (r, &v) in unit.iter().enumerate() {
                    inverse[r * nr + c] = v;
                }
            }
        }
        // Back substitution: x_j = S_j^{-1} (z_j + A_j x_{j+1}).
        let mut response = vec![0.0; nw * nw];
        let mut next = vec![0.0; nr * nw];
        let mut work = vec![0.0; nr * nw];
        for j in (0..nx).rev() {
            work.copy_from_slice(&z[j * nr * nw..(j + 1) * nr * nw]);
            if j + 1 < nx {
                let ax = &self.axial[j * nr..(j + 1) * nr];
                for (p, &g) in ax.iter().enumerate() {
                    for (o, &v) in work[p * nw..(p + 1) * nw]
                        .iter_mut()
                        .zip(&next[p * nw..(p + 1) * nw])
                    {
                        *o += g * v;
                    }
                }
            }
            let inverse = &inverses[j * block..(j + 1) * block];
            next.fill(0.0);
            for p in 0..nr {
                let out = &mut next[p * nw..(p + 1) * nw];
                for q in 0..nr {
                    let f = inverse[p * nr + q];
                    for (o, &v) in out.iter_mut().zip(&work[q * nw..(q + 1) * nw]) {
                        *o += f * v;
                    }
                }
            }
            if self.columns.contains(&j) {
                response[(j - c0) * nw..(j - c0 + 1) * nw]
                    .copy_from_slice(&next[row * nw..(row + 1) * nw]);
            }
        }
        Ok(response)
    }
}

/// Edges from `x0` to `x1` (either order): cells of `h0` at `x0`, growing by
/// `ratio`, the last one absorbing the remainder.
fn graded(x0: f64, x1: f64, h0: f64, ratio: f64, cap: usize) -> Result<Vec<f64>> {
    let span = (x1 - x0).abs();
    let sign = if x1 >= x0 { 1.0 } else { -1.0 };
    let mut widths = Vec::new();
    let (mut covered, mut h) = (0.0, h0);
    while covered + h < span {
        if widths.len() == cap {
            return Err(WallWindowError::Budget {
                what: "cells along a grid line",
                required: cap + 1,
                maximum: cap,
            });
        }
        widths.push(h);
        covered += h;
        h *= ratio;
    }
    match widths.last_mut() {
        Some(last) => *last += span - covered,
        None => widths.push(span),
    }
    let mut edges = Vec::with_capacity(widths.len() + 1);
    let mut x = x0;
    edges.push(x);
    for width in widths {
        x += sign * width;
        edges.push(x);
    }
    let n = edges.len();
    edges[n - 1] = x1;
    Ok(edges)
}

/// Orthonormal DCT-II rows, row-major `n x n`.
fn dct_rows(n: usize) -> Vec<f64> {
    let mut rows = vec![0.0; n * n];
    for p in 0..n {
        let norm = det::sqrt(if p == 0 { 1.0 } else { 2.0 } / n as f64);
        for k in 0..n {
            rows[p * n + k] = norm * det::cos(PI * p as f64 * (k as f64 + 0.5) / n as f64);
        }
    }
    rows
}

/// `P Y P^T` for row-major `n x n` matrices.
fn congruence(p: &[f64], y: &[f64], n: usize) -> Vec<f64> {
    let mut py = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let f = p[i * n + k];
            for (o, &v) in py[i * n..(i + 1) * n]
                .iter_mut()
                .zip(&y[k * n..(k + 1) * n])
            {
                *o += f * v;
            }
        }
    }
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            out[i * n + j] = py[i * n..(i + 1) * n]
                .iter()
                .zip(&p[j * n..(j + 1) * n])
                .map(|(x, y)| x * y)
                .sum();
        }
    }
    out
}

/// Corners `(K^-1)_00`, `(K^-1)_0L` and `(K^-1)_LL` of the channel's
/// negated radial operator in one DCT mode. It is tridiagonal over the wall
/// layers, with the two face potentials half a layer beyond the end layers.
fn channel_inverse_corners(
    edge_r: &[f64],
    layer: f64,
    mu_p: f64,
    mu_q: f64,
    gb: f64,
    gt: f64,
) -> (f64, f64, f64) {
    let n = edge_r.len() - 1;
    let off: Vec<f64> = (1..n).map(|l| edge_r[l] / layer).collect();
    let diag: Vec<f64> = (0..n)
        .map(|l| {
            let rc = 0.5 * (edge_r[l] + edge_r[l + 1]);
            let mut d = mu_p * layer / rc + mu_q * rc * layer;
            if l > 0 {
                d += off[l - 1];
            }
            if l + 1 < n {
                d += off[l];
            }
            if l == 0 {
                d += gb;
            }
            if l + 1 == n {
                d += gt;
            }
            d
        })
        .collect();
    let solve = |unit: usize| -> Vec<f64> {
        let mut c = vec![0.0; n];
        let mut x = vec![0.0; n];
        let mut pivot = diag[0];
        x[0] = if unit == 0 { 1.0 } else { 0.0 } / pivot;
        for l in 1..n {
            c[l - 1] = -off[l - 1] / pivot;
            pivot = diag[l] + off[l - 1] * c[l - 1];
            let b = if unit == l { 1.0 } else { 0.0 };
            x[l] = (b + off[l - 1] * x[l - 1]) / pivot;
        }
        for l in (0..n - 1).rev() {
            x[l] -= c[l] * x[l + 1];
        }
        x
    };
    let first = solve(0);
    let last = solve(n - 1);
    (first[0], first[n - 1], last[n - 1])
}
