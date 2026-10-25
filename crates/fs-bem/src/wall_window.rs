//! Potential-flow acoustic mass of a window cut through the wall of a
//! cylindrical bore: a flue instrument's mouth as the air column sees it.
//!
//! The bore has radius `a` and a rigid wall of thickness `t`. The window is
//! cut radially through the wall over the arc `|theta| <= theta_w` and the
//! axial interval `|x| <= L/2`, and its side walls are rigid. Upstream, the
//! bore ends in a flat closed end (a flue's pocket) or in a port; downstream
//! it continues to a port. Outside, the wall is an infinite rigid cylinder of
//! radius `a + t` in free space. Alternatively the window's outer face is a
//! pressure release, which leaves no outer region and is used to check the
//! inner end correction on its own.
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
//! - The channel has `wall_layers` radial layers and is diagonal in the
//!   window cells' DCT-II basis. Per mode, it maps its two face potentials to
//!   their fluxes.
//! - The three regions couple through the two faces in that basis.
//!   Geometry and forcing are symmetric about `theta = 0`, so only the
//!   `theta`-even DCT modes enter, and each bore mode's contribution has rank
//!   one in `theta`. The coupled system has `2 ceil(arc_cells / 2)
//!   axial_cells` unknowns (half that with a release face) and is solved by
//!   LU.
//! - The exterior's far boundary is a monopole Robin condition about the
//!   origin. In the cases checked, its error falls about eightfold for each
//!   doubling of the far radius.
//!
//! The arc is resolved as `arc_cells` whole cells of `2 pi / N`, with `N` the
//! nearest integer. The resolved arc can therefore differ from the requested
//! one; it is reported, and a ladder should use rungs that resolve it alike.
//! The window's edges are re-entrant right angles, where the velocity is
//! singular like `d^(-1/3)`, so the mass converges at order 4/3 in the window
//! cell size (1.30 observed). A caller ladders the cells and extrapolates.
//! Nothing is fitted.
//!
//! Not represented: window walls that are not radial; a flat baffle,
//! windway or labium around the window; a pocket that is not flat;
//! viscosity; compressibility.

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

/// A window cut radially through the wall of a cylindrical bore.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallWindow {
    /// Bore radius `a` [m].
    pub bore_radius: f64,
    /// Wall thickness `t` at the window [m].
    pub wall_thickness: f64,
    /// Window width as arc length on the bore surface, `2 theta_w a` [m].
    pub arc_width: f64,
    /// Window length along the bore axis [m].
    pub axial_length: f64,
    /// Upstream closure.
    pub upstream: WindowUpstream,
    /// Outer side.
    pub outside: WindowOutside,
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
    /// A coupled face system with more unknowns than this is refused.
    pub max_unknowns: usize,
    /// Workers for the `theta` modes. Results are bitwise identical at every
    /// count.
    pub workers: usize,
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
        }
    }
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
    /// Resolved window area `S_w` [m^2].
    pub window_area: f64,
    /// Resolved arc width [m].
    pub resolved_arc_width: f64,
    /// `Q`: flux into the window's bore face per unit potential [m].
    pub flux: f64,
    /// Flux out of the window's outer face per unit potential [m]. It equals
    /// `flux` to rounding.
    pub outer_flux: f64,
    /// `theta` cells around the bore, `N`.
    pub theta_cells: usize,
    /// Bore grid cells, (radial, axial).
    pub bore_cells: (usize, usize),
    /// Exterior grid cells, (radial, axial); zero with a release face.
    pub exterior_cells: (usize, usize),
    /// Unknowns of the coupled face system.
    pub unknowns: usize,
    /// Relative residual of the coupled face system, recomputed after the
    /// solve.
    pub residual: f64,
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
    let even = nwth.div_ceil(2);
    let face = even * nwx;
    let unknowns = if free { 2 * face } else { face };
    if unknowns > o.max_unknowns {
        return Err(WallWindowError::Budget {
            what: "coupled face unknowns",
            required: unknowns,
            maximum: o.max_unknowns,
        });
    }

    let dth = 2.0 * PI / n_theta as f64;
    let dxw = len / nwx as f64;
    let h = dxw.min(a * dth);
    let (ratio, cap) = (o.grading, o.max_line_cells);
    let window_edges: Vec<f64> = (0..=nwx)
        .map(|i| -half + len * i as f64 / nwx as f64)
        .collect();
    let axial = |far_up: f64, far_down: f64, pocket: bool| -> Result<(Vec<f64>, usize)> {
        let mut edges = Vec::new();
        if pocket {
            let mut up = graded(-half, -far_up, dxw, ratio, cap)?;
            up.reverse();
            edges.extend_from_slice(&up[..up.len() - 1]);
        }
        let first = edges.len();
        edges.extend_from_slice(&window_edges);
        edges.extend_from_slice(&graded(half, far_down, dxw, ratio, cap)?[1..]);
        Ok((edges, first))
    };
    let (bore_x, bore_first) = axial(upstream, port, upstream - half > 1e-6 * dxw)?;
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
        let (ext_x, ext_first) = axial(reach, reach, true)?;
        let ext_r = graded(a + t, reach, h, ratio, cap)?;
        Some(Region::new(
            &ext_r,
            &ext_x,
            0,
            ext_first..ext_first + nwx,
            Far::Robin,
        ))
    } else {
        None
    };

    let modes = n_theta / 2 + 1;
    let context = Modes {
        n_theta,
        dth,
        arc_cells: nwth,
        axial_cells: nwx,
        even,
        arc_dct: dct_rows(nwth),
        axial_dct: dct_rows(nwx),
        bore: &bore,
        exterior: exterior.as_ref(),
    };
    let per = modes.div_ceil(o.workers.min(modes));
    let chunks: Vec<Result<Vec<ModeTerms>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..modes)
            .step_by(per)
            .map(|start| {
                let context = &context;
                scope.spawn(move || {
                    (start..(start + per).min(modes))
                        .map(|m| context.terms(m))
                        .collect()
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("wall window mode worker panicked"))
            .collect()
    });
    let mut terms = Vec::with_capacity(modes);
    for chunk in chunks {
        terms.extend(chunk?);
    }

    // Face maps in the even DCT basis: potential per unit flux out of the
    // region through each mode, plus the face cells' half-cell drop.
    let face_map = |pick: fn(&ModeTerms) -> &[f64], half_cell: f64, radius: f64| {
        let mut map = vec![0.0; face * face];
        for (m, term) in terms.iter().enumerate() {
            let weight = if m == 0 || 2 * m == n_theta { 1.0 } else { 2.0 };
            let g = weight / (n_theta as f64 * dth);
            let y = pick(term);
            for pe in 0..even {
                for pf in 0..even {
                    let coefficient = g * term.arc[pe] * term.arc[pf];
                    for q in 0..nwx {
                        let row = &mut map[(pe * nwx + q) * face + pf * nwx..][..nwx];
                        for (entry, &yv) in row.iter_mut().zip(&y[q * nwx..(q + 1) * nwx]) {
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
    let bore_map = face_map(bore_terms, bore.half_cell, a);
    let exterior_map = exterior
        .as_ref()
        .map(|region| face_map(exterior_terms, region.half_cell, a + t));

    // Channel: per even DCT mode, flux into the bottom face and out of the
    // top face per unit bottom and top potential, scaled to cell fluxes.
    let mut into_bottom = vec![0.0; face];
    let mut cross = vec![0.0; face];
    let mut out_top = vec![0.0; face];
    let layer = t / o.wall_layers as f64;
    let edge_r: Vec<f64> = (0..=o.wall_layers)
        .map(|l| a + t * l as f64 / o.wall_layers as f64)
        .collect();
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
    let size = unknowns;
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
    let residual = det::sqrt(misfit) / rhs[0];

    let top = if free { phi[face] } else { 0.0 };
    let flux = rhs[0] * (into_bottom[0] * phi[0] + cross[0] * top);
    let outer_flux = rhs[0] * (-cross[0] * phi[0] + out_top[0] * top);
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
        resolved_arc_width,
        flux,
        outer_flux,
        theta_cells: n_theta,
        bore_cells: (bore.nr, bore.nx),
        exterior_cells: exterior.as_ref().map_or((0, 0), |e| (e.nr, e.nx)),
        unknowns,
        residual,
    })
}

/// What every `theta` mode's terms share.
struct Modes<'a> {
    n_theta: usize,
    dth: f64,
    arc_cells: usize,
    axial_cells: usize,
    even: usize,
    /// Orthonormal DCT-II rows over the window's arc cells.
    arc_dct: Vec<f64>,
    /// Orthonormal DCT-II rows over the window's axial cells.
    axial_dct: Vec<f64>,
    bore: &'a Region,
    exterior: Option<&'a Region>,
}

impl Modes<'_> {
    /// Mode `m`: the discrete `theta` Laplacian's eigenvalue `mu`, the arc
    /// coefficients of `cos(m theta)` sampled at the arc cells' centres, and
    /// each region's face response projected onto the axial DCT basis.
    fn terms(&self, m: usize) -> Result<ModeTerms> {
        let (nwth, nwx, n) = (self.arc_cells, self.axial_cells, self.n_theta as f64);
        let mu = (2.0 - 2.0 * det::cos(2.0 * PI * m as f64 / n)) / (self.dth * self.dth);
        let cos_m: Vec<f64> = (0..nwth)
            .map(|k| {
                let centre = k as f64 - 0.5 * (nwth - 1) as f64;
                det::cos(2.0 * PI * m as f64 * centre / n)
            })
            .collect();
        let arc = (0..self.even)
            .map(|pe| {
                let row = &self.arc_dct[2 * pe * nwth..(2 * pe + 1) * nwth];
                row.iter().zip(&cos_m).map(|(p, c)| p * c).sum::<f64>()
            })
            .collect();
        let response = |region: &Region, stage: &'static str| {
            region
                .face_response(mu)
                .map(|y| congruence(&self.axial_dct, &y, nwx))
                .map_err(|source| WallWindowError::Factor { stage, source })
        };
        let bore = response(self.bore, "bore mode elimination")?;
        let exterior = match self.exterior {
            Some(region) => response(region, "exterior mode elimination")?,
            None => Vec::new(),
        };
        Ok(ModeTerms {
            arc,
            bore,
            exterior,
        })
    }
}

/// One `theta` mode's contribution: the even arc DCT coefficients of
/// `cos(m theta)` over the window's cells, and each region's face response
/// in the axial DCT basis (`axial_cells` squared, row-major).
struct ModeTerms {
    arc: Vec<f64>,
    bore: Vec<f64>,
    exterior: Vec<f64>,
}

fn bore_terms(term: &ModeTerms) -> &[f64] {
    &term.bore
}

fn exterior_terms(term: &ModeTerms) -> &[f64] {
    &term.exterior
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
