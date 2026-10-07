//! `frankensim cooling-cht <scene.json>`: file-driven voxel conjugate heat
//! transfer. An explicit experimental JSON workflow (like `cooling-network`),
//! not a `.fsim` or ledger-backed solve: the scene declares a box domain at
//! one voxel size, solid boxes with conductivities, heat-source boxes, and
//! one flow/thermal rule per domain face; the command solves steady laminar
//! flow by finite-volume SIMPLEC (forced convection) or SIMPLEC with the
//! Boussinesq force (natural/mixed convection, when `gravity_m_s2` is
//! declared) and the conservative conjugate energy equation over fluid and
//! solid cells (fs-lbm `conjugate`).

#[path = "json_read.rs"]
#[allow(dead_code)]
mod json;
mod study;

use std::ffi::OsString;
use std::fmt::{self, Write as _};
use std::fs::File;
use std::io::Read;
use std::time::{Duration, Instant};

use fs_cli::{CommandOutput, exit};
use fs_exec::CancelGate;
use fs_geom::Point3;
use fs_io::stl::read_stl;
use fs_lbm::Face3;
use fs_lbm::conjugate::{
    Boussinesq, ChtError, CompactComponent, ContactResistance, EnergyConfig, EnergySolution,
    FacePatch, FanCurve, FanInlet, FlowResistance, FluidProperties, FvBoundary, FvBuoyancyConfig,
    FvFlow, InternalFan, RadiationConfig, SimpleConfig, SolidMaterial, ThermalFace, ThermalSetup,
    TimeScheme, TransientConfig, Turbulence, UnsteadyConfig, Voxel, VoxelDomain,
    fv_natural_convection, march_conjugate, march_energy, simple_flow, simple_unsteady,
    solve_energy, solve_energy_radiating,
};
use fs_rep_mesh::{Soup, WindingOctree, winding_exact};
use json::JsonValue as J;

const MAX_INPUT_BYTES: u64 = 16 * 1024 * 1024;
const MAX_CELLS: usize = 4_000_000;
const MAX_STL_BYTES: u64 = 256 * 1024 * 1024;
const SCHEMA: &str = "frankensim.cooling-cht.v1";
const RESULT_SCHEMA: &str = "frankensim.cooling-cht.result.v1";
const NO_CLAIM: &str = "steady constant-property flow on a staircase voxel grid at one declared resolution (no mesh-convergence claim); Boussinesq buoyancy only when gravity is declared; turbulence only through the optional LVEL algebraic eddy viscosity (no transport, separation or transition physics); no temperature-dependent properties; radiation only between gray diffuse exposed solid faces and to the surroundings seen through openings, inlets and fans (Monte Carlo exchange factors on face patches; walls and non-emitting solids reflect perfectly; transparent air); power-law convection is first order at high cell Peclet numbers; Estimated numerical evidence, not validated hardware or a ledger-backed .fsim run";
const HELP: &str = "Usage: frankensim [--json] cooling-cht <scene.json>\n\nSolve steady voxel conjugate heat transfer: finite-volume SIMPLEC airflow\n(forced, or natural/mixed with the Boussinesq force when gravity_m_s2 is\ndeclared) and one conservative energy equation over fluid and solid cells.\nThe scene declares size_m and voxel_m, a fluid (\"dry-air-300k\" or explicit\nproperties), materials (isotropic k or [kx, ky, kz]), contacts (interface\nresistance_m2_k_w between two materials), solids (boxes, or closed STL meshes placed by\nscale and offset_m; later solids override earlier ones),\nheat-source boxes (power spread over the solid cells they cover), and one\nrule per face x-, x+, y-, y+, z-, z+: inlet (velocity_m_s, temperature_k),\nfan (curve [[flow_m3_s, pressure_pa], ...], temperature_k; the flow is the\noperating point against the system), opening (ambient_k; pressure zero, flow either way), symmetry, or wall\n(adiabatic, or temperature_k, heat_flux_w_m2, or htc_w_m2_k with ambient_k).\nMissing faces are adiabatic walls. A material emissivity enables gray\nsurface radiation between emitting faces and to the surroundings seen\nthrough openings, inlets and fans (Monte Carlo exchange factors; walls\nand non-emitting solids reflect; radiation {rays_per_face, seed,\nsurface_exchange (default true), patch_size (default 4)}).\nsolver.turbulence \"lvel\" adds the LVEL algebraic eddy viscosity (and its\nturbulent conductivity) for transitional/turbulent fan-driven flow.\ninternal_fans (axis, at_m on an interior voxel face, direction \"+\"/\"-\",\nmin_m/max_m transverse extent, curve) raise the pressure across a plane;\nresistances are grilles (axis, at_m, min_m/max_m, loss_coefficient or\nfree_area_ratio) or porous blocks (min_m/max_m, permeability_m2 and\ninertial_per_m, scalar or per axis).\nA solid may be a plate-fin heatsink (heatsink {base_min_m, base_size_m,\nfin_count, fin_thickness_m, fin_height_m, fins_along}). A study block\n(parameters [{name, path, values}], objective {minimize}, constraints\n[{quantity, min, max}]) evaluates every combination of the values and ranks\nthe variants; quantities are max_solid_temperature_k, source:<name>,\ncomponent:<name>, internal_fan:<name>, fan_flow_m3_s and inflow_m3_s.\ncomponents are JEDEC two-resistor compact models (min_m/max_m box,\nboard_side, power_w, junction_to_case_k_w, junction_to_board_k_w): the box\nblocks flow and the junction reaches the case top and the board through\nthe two resistors (steady scenes only).\nOptional transient (time_step_s, steps,\npower_schedule [[time_s, scale], ...], initial_temperature_k) marches the\nenergy equation over the steady forced flow (materials then need\nvolumetric_heat_capacity_j_m3_k); with flow \"unsteady\" (scheme \"bdf2\" or\n\"backward-euler\", inner_iterations, inner_tolerance, inlet_schedule) the\nflow marches with it from rest, buoyant when gravity is declared; energy\n\"steady-on-mean-flow\" instead solves the steady energy equation (with any\nradiation) on the march's time-averaged fluxes. Request schema: frankensim.cooling-cht.v1.\nResults are Estimated single-resolution numerical evidence.\n";

type Result<T> = std::result::Result<T, Failure>;

#[derive(Debug)]
struct Failure {
    code: &'static str,
    message: String,
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}

fn bad(message: impl Into<String>) -> Failure {
    Failure {
        code: "cooling-cht-input",
        message: message.into(),
    }
}

fn solver_failure(error: &ChtError) -> Failure {
    let code = match error {
        ChtError::Cancelled => "cooling-cht-cancelled",
        ChtError::FlowNotSteady { .. } => "cooling-cht-not-steady",
        ChtError::InvalidInput { .. } | ChtError::InvalidDomain { .. } => "cooling-cht-input",
        _ => "cooling-cht-solve",
    };
    Failure {
        code,
        message: error.to_string(),
    }
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c < '\u{20}' => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// A finite number, or `null` (a quantity the scene does not have, such as
/// a solid temperature without solids).
fn num_or_null(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "null".into()
    }
}

fn num(value: f64) -> Result<String> {
    if value.is_finite() {
        Ok(value.to_string())
    } else {
        Err(Failure {
            code: "cooling-cht-solve",
            message: "nonfinite result cannot be published as JSON".into(),
        })
    }
}

fn field<'a>(object: &'a J, key: &str, at: &str) -> Result<&'a J> {
    object
        .get(key)
        .ok_or_else(|| bad(format!("{at}.{key} is required")))
}

fn number(object: &J, key: &str, at: &str) -> Result<f64> {
    let value = field(object, key, at)?
        .as_f64()
        .ok_or_else(|| bad(format!("{at}.{key} must be a number")))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(format!("{at}.{key} must be finite")))
    }
}

fn optional_number(object: &J, key: &str, at: &str) -> Result<Option<f64>> {
    match object.get(key) {
        None => Ok(None),
        Some(_) => number(object, key, at).map(Some),
    }
}

fn vec3(object: &J, key: &str, at: &str) -> Result<[f64; 3]> {
    let items = field(object, key, at)?
        .as_array()
        .filter(|items| items.len() == 3)
        .ok_or_else(|| bad(format!("{at}.{key} must be a 3-vector")))?;
    let mut out = [0.0; 3];
    for (slot, item) in out.iter_mut().zip(items) {
        *slot = item
            .as_f64()
            .filter(|v| v.is_finite())
            .ok_or_else(|| bad(format!("{at}.{key} entries must be finite numbers")))?;
    }
    Ok(out)
}

/// The array under `key` (absent: empty).
fn array_of<'a>(object: &'a J, key: &str) -> Result<&'a [J]> {
    match object.get(key) {
        None => Ok(&[]),
        Some(value) => value
            .as_array()
            .ok_or_else(|| bad(format!("{key} must be an array"))),
    }
}

/// A fan characteristic `curve: [[flow_m3_s, pressure_pa], ...]`.
fn fan_curve(rule: &J, at: &str) -> Result<FanCurve> {
    let points = field(rule, "curve", at)?
        .as_array()
        .ok_or_else(|| {
            bad(format!(
                "{at}.curve must be [[flow_m3_s, pressure_pa], ...]"
            ))
        })?
        .iter()
        .map(|point| {
            point
                .as_array()
                .filter(|pair| pair.len() == 2)
                .and_then(|pair| Some((pair[0].as_f64()?, pair[1].as_f64()?)))
                .ok_or_else(|| {
                    bad(format!(
                        "{at}.curve entries must be [flow_m3_s, pressure_pa]"
                    ))
                })
        })
        .collect::<Result<Vec<_>>>()?;
    FanCurve::new(&points).map_err(|e| bad(format!("{at}.curve: {e}")))
}

/// A number or a 3-vector under `key` (absent: `default` on every axis).
fn per_axis(object: &J, key: &str, at: &str, default: f64) -> Result<[f64; 3]> {
    match object.get(key) {
        None => Ok([default; 3]),
        Some(value) if value.as_f64().is_some() => Ok([number(object, key, at)?; 3]),
        Some(_) => vec3(object, key, at),
    }
}

/// The scene's Cartesian grid: uniform (`voxel_m`) or graded per axis
/// (`grid` zones).
#[derive(Debug, Clone)]
struct Grid {
    widths: [Vec<f64>; 3],
    faces: [Vec<f64>; 3],
    centres: [Vec<f64>; 3],
    /// The uniform spacing, when uniform.
    uniform: Option<f64>,
    min_width: f64,
}

impl Grid {
    fn from_widths(widths: [Vec<f64>; 3]) -> Self {
        let faces = widths.clone().map(|list| {
            let mut faces = vec![0.0];
            for w in list {
                faces.push(faces.last().copied().unwrap_or(0.0) + w);
            }
            faces
        });
        let centres = faces
            .clone()
            .map(|f| f.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect());
        let first = widths[0][0];
        let uniform = widths
            .iter()
            .all(|list| list.iter().all(|w| *w == first))
            .then_some(first);
        let min_width = widths
            .iter()
            .flatten()
            .fold(f64::INFINITY, |m, w| m.min(*w));
        Self {
            widths,
            faces,
            centres,
            uniform,
            min_width,
        }
    }

    fn dims(&self) -> [usize; 3] {
        [0, 1, 2].map(|a| self.widths[a].len())
    }

    fn length(&self, axis: usize) -> f64 {
        *self.faces[axis].last().expect("non-empty axis")
    }

    /// Cells of `axis` whose tie-shifted centres lie in `[min, max)`.
    fn span(&self, axis: usize, min: f64, max: f64) -> (usize, usize) {
        let tie = 1e-6 * self.min_width;
        let centres = &self.centres[axis];
        let first = |edge: f64| centres.partition_point(|c| c + tie < edge);
        (first(min), first(max))
    }

    /// The face plane of `axis` at coordinate `at` (within 1e-6 of the
    /// smallest width).
    fn plane(&self, axis: usize, at: f64) -> Option<usize> {
        let faces = &self.faces[axis];
        let i = faces.partition_point(|f| *f < at - 1e-6 * self.min_width);
        (i < faces.len() && (faces[i] - at).abs() <= 1e-6 * self.min_width).then_some(i)
    }

    /// The cell of `axis` containing coordinate `v`.
    fn locate(&self, axis: usize, v: f64) -> usize {
        self.faces[axis]
            .partition_point(|f| *f <= v)
            .clamp(1, self.widths[axis].len())
            - 1
    }

    fn domain(
        &self,
        occupancy: impl FnMut([f64; 3]) -> Voxel,
    ) -> std::result::Result<VoxelDomain, ChtError> {
        match self.uniform {
            Some(dx) => {
                let [nx, ny, nz] = self.dims();
                VoxelDomain::from_fn(nx, ny, nz, dx, occupancy)
            }
            None => VoxelDomain::graded_from_fn(self.widths.clone(), occupancy),
        }
    }
}

/// Cell ranges `lo..hi` of a box (each non-empty) in the grid.
fn cell_box(region: &Aabb, grid: &Grid, at: &str) -> Result<([usize; 3], [usize; 3])> {
    let (mut lo, mut hi) = ([0; 3], [0; 3]);
    for a in 0..3 {
        (lo[a], hi[a]) = grid.span(a, region.min[a], region.max[a]);
        if lo[a] >= hi[a] {
            return Err(bad(format!(
                "{at}: the box covers no voxel centre on axis {a}"
            )));
        }
    }
    Ok((lo, hi))
}

/// An interior face patch: `axis` ("x", "y", "z"), the plane `at_m` (on a
/// voxel face strictly inside the domain), and the transverse extent of the
/// box `min_m`/`max_m` (its entries along `axis` are ignored).
fn face_patch(item: &J, grid: &Grid, at: &str) -> Result<FacePatch> {
    let dims = grid.dims();
    let axis = match item.str_field("axis") {
        Some("x") => 0,
        Some("y") => 1,
        Some("z") => 2,
        _ => return Err(bad(format!("{at}.axis must be \"x\", \"y\" or \"z\""))),
    };
    let index = grid
        .plane(axis, number(item, "at_m", at)?)
        .filter(|&i| i >= 1 && i < dims[axis])
        .ok_or_else(|| {
            bad(format!(
                "{at}.at_m must lie on a voxel face strictly inside the domain"
            ))
        })?;
    let mut min = vec3(item, "min_m", at)?;
    let mut max = vec3(item, "max_m", at)?;
    // The normal extent is irrelevant: give the box one full cell there.
    min[axis] = 0.0;
    max[axis] = grid.faces[axis][1];
    let (lo, hi) = cell_box(&Aabb { min, max }, grid, at)?;
    Ok(FacePatch {
        axis,
        index,
        lo,
        hi,
    })
}

/// A plate-fin heatsink: a base box `base_min_m` + `base_size_m` and
/// `fin_count` plates of `fin_thickness_m` standing `fin_height_m` on its
/// top (+z), running along `fins_along` (`"x"` or `"y"`) and spread evenly
/// across the other horizontal axis (the outer fins flush with the base
/// edges). Every fin must cover voxel centres and every gap must keep a
/// fluid voxel, so a design the grid cannot represent refuses.
fn plate_fin_heatsink(item: &J, grid: &Grid, at: &str) -> Result<Vec<Aabb>> {
    let origin = vec3(item, "base_min_m", at)?;
    let size = vec3(item, "base_size_m", at)?;
    let count = number(item, "fin_count", at)?;
    let thickness = number(item, "fin_thickness_m", at)?;
    let height = number(item, "fin_height_m", at)?;
    if size.iter().any(|v| *v <= 0.0)
        || !(1.0..=1000.0).contains(&count)
        || count.fract() != 0.0
        || thickness <= 0.0
        || height <= 0.0
    {
        return Err(bad(format!(
            "{at}: base_size_m, fin_thickness_m and fin_height_m must be positive and fin_count a whole number in 1..=1000"
        )));
    }
    // The axis the fins are spread across (they run along the other).
    let across = match item.str_field("fins_along") {
        None | Some("x") => 1,
        Some("y") => 0,
        Some(_) => return Err(bad(format!("{at}.fins_along must be \"x\" or \"y\""))),
    };
    let count = count as usize;
    if thickness * count as f64 > size[across] + 1e-12 {
        return Err(bad(format!(
            "{at}: {count} fins of {thickness} m do not fit the base"
        )));
    }
    let mut parts = vec![Aabb {
        min: origin,
        max: [0, 1, 2].map(|a| origin[a] + size[a]),
    }];
    let pitch = if count > 1 {
        (size[across] - thickness) / (count - 1) as f64
    } else {
        0.0
    };
    let mut previous_end: Option<f64> = None;
    for i in 0..count {
        let start = if count > 1 {
            origin[across] + i as f64 * pitch
        } else {
            origin[across] + 0.5 * (size[across] - thickness)
        };
        let mut min = origin;
        let mut max = [0, 1, 2].map(|a| origin[a] + size[a]);
        min[across] = start;
        max[across] = start + thickness;
        min[2] = origin[2] + size[2];
        max[2] = min[2] + height;
        let (lo, hi) = grid.span(across, min[across], max[across]);
        if lo >= hi {
            return Err(bad(format!(
                "{at}: fin {i} ({thickness} m) covers no voxel centre on this grid"
            )));
        }
        if let Some(end) = previous_end {
            let (gap_lo, gap_hi) = grid.span(across, end, start);
            if gap_lo >= gap_hi {
                return Err(bad(format!(
                    "{at}: the gap before fin {i} holds no fluid voxel on this grid (fins merge)"
                )));
            }
        }
        previous_end = Some(start + thickness);
        parts.push(Aabb { min, max });
    }
    Ok(parts)
}

/// A graded grid: per axis (`x`, `y`, `z`) a list of zones `{"to_m",
/// "voxel_m"}` from the previous zone's end (0 first), each a whole number
/// of uniform cells.
fn parse_grid(value: &J) -> Result<Grid> {
    let mut widths = [Vec::new(), Vec::new(), Vec::new()];
    for (a, key) in ["x", "y", "z"].iter().enumerate() {
        let at = format!("grid.{key}");
        let zones = value
            .get(key)
            .and_then(J::as_array)
            .filter(|z| !z.is_empty())
            .ok_or_else(|| {
                bad(format!(
                    "{at} must be a non-empty array of {{to_m, voxel_m}} zones"
                ))
            })?;
        let mut from = 0.0f64;
        for (i, zone) in zones.iter().enumerate() {
            let zat = format!("{at}[{i}]");
            let to = number(zone, "to_m", &zat)?;
            let voxel = number(zone, "voxel_m", &zat)?;
            if !(to > from) || !(voxel > 0.0) {
                return Err(bad(format!(
                    "{zat}: to_m must increase and voxel_m be positive"
                )));
            }
            let cells = (to - from) / voxel;
            let rounded = cells.round();
            if rounded < 1.0 || (cells - rounded).abs() > 1e-6 * rounded.max(1.0) {
                return Err(bad(format!(
                    "{zat}: the zone {from}..{to} m is not a whole number of {voxel} m voxels"
                )));
            }
            let width = (to - from) / rounded;
            widths[a].extend(std::iter::repeat_n(width, rounded as usize));
            from = to;
        }
    }
    if widths.iter().map(Vec::len).product::<usize>() > MAX_CELLS {
        return Err(Failure {
            code: "cooling-cht-budget",
            message: format!(
                "{} cells exceed the {MAX_CELLS}-cell cap",
                widths.iter().map(Vec::len).product::<usize>()
            ),
        });
    }
    Ok(Grid::from_widths(widths))
}

/// An axis-aligned box `[min, max)` in metres.
#[derive(Debug, Clone, Copy)]
struct Aabb {
    min: [f64; 3],
    max: [f64; 3],
}

impl Aabb {
    fn parse(object: &J, at: &str) -> Result<Self> {
        let min = vec3(object, "min_m", at)?;
        let max = vec3(object, "max_m", at)?;
        if (0..3).any(|a| min[a] >= max[a]) {
            return Err(bad(format!(
                "{at}: min_m must be below max_m on every axis"
            )));
        }
        Ok(Self { min, max })
    }

    /// Membership of a query point (see [`probe`]).
    fn contains(&self, q: [f64; 3]) -> bool {
        (0..3).all(|a| q[a] >= self.min[a] && q[a] < self.max[a])
    }
}

/// The point at which a voxel centre's membership is evaluated: the centre
/// moved up by `tie` (a tiny fraction of the voxel) on every axis, so a box
/// edge or mesh face lying exactly on a centre deterministically takes the
/// voxel at a box's min edge and leaves the one at its max edge, whatever
/// the rounding of decimal coordinates.
fn probe(centre: [f64; 3], tie: f64) -> [f64; 3] {
    centre.map(|v| v + tie)
}

/// A closed triangle mesh placed by `world = scale * stl + offset`; inside
/// is the robust generalized winding number above one half.
struct MeshSolid {
    soup: Soup,
    tree: Option<WindingOctree>,
    scale: f64,
    offset: [f64; 3],
    bounds: Aabb,
}

/// Triangle count above which classification uses the dipole octree
/// instead of the exact solid-angle sum.
const EXACT_WINDING_TRIANGLES: usize = 4096;

impl MeshSolid {
    fn load(item: &J, at: &str, base: &std::path::Path) -> Result<Self> {
        let path = base.join(
            item.str_field("stl")
                .ok_or_else(|| bad(format!("{at}.stl must be a path string")))?,
        );
        let scale = optional_number(item, "scale", at)?.unwrap_or(1.0);
        if scale <= 0.0 {
            return Err(bad(format!("{at}.scale must be positive")));
        }
        let offset = match item.get("offset_m") {
            None => [0.0; 3],
            Some(_) => vec3(item, "offset_m", at)?,
        };
        let mut bytes = Vec::new();
        File::open(&path)
            .and_then(|file| file.take(MAX_STL_BYTES + 1).read_to_end(&mut bytes))
            .map_err(|e| bad(format!("{at}: cannot read {}: {e}", path.display())))?;
        if bytes.len() as u64 > MAX_STL_BYTES {
            return Err(bad(format!("{at}: STL exceeds {MAX_STL_BYTES} bytes")));
        }
        let soup = read_stl(&bytes).map_err(|e| bad(format!("{at}: {}: {e:?}", path.display())))?;
        if soup.triangles.is_empty() {
            return Err(bad(format!("{at}: STL has no triangles")));
        }
        let mut bounds = Aabb {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        };
        for p in &soup.positions {
            for (a, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                bounds.min[a] = bounds.min[a].min(v);
                bounds.max[a] = bounds.max[a].max(v);
            }
        }
        let tree = (soup.triangles.len() > EXACT_WINDING_TRIANGLES)
            .then(|| WindingOctree::build(&soup, 2.0));
        Ok(Self {
            soup,
            tree,
            scale,
            offset,
            bounds,
        })
    }

    fn contains(&self, q: [f64; 3]) -> bool {
        let local = [0, 1, 2].map(|a| (q[a] - self.offset[a]) / self.scale);
        if (0..3).any(|a| local[a] < self.bounds.min[a] || local[a] > self.bounds.max[a]) {
            return false;
        }
        let point = Point3 {
            x: local[0],
            y: local[1],
            z: local[2],
        };
        match &self.tree {
            Some(tree) => tree.inside(&self.soup, point),
            None => winding_exact(&self.soup, point) > 0.5,
        }
    }
}

enum Shape {
    Box(Aabb),
    Mesh(Box<MeshSolid>),
}

impl Shape {
    fn contains(&self, q: [f64; 3]) -> bool {
        match self {
            Self::Box(region) => region.contains(q),
            Self::Mesh(mesh) => mesh.contains(q),
        }
    }
}

#[derive(Debug)]
struct Source {
    name: String,
    power_w: f64,
    region: Aabb,
}

#[derive(Debug, Clone, Copy)]
enum FaceRule {
    Inlet {
        velocity: [f64; 3],
        temperature: f64,
    },
    Fan {
        curve: FanCurve,
        temperature: f64,
    },
    Opening {
        ambient: f64,
    },
    Symmetry,
    Wall(ThermalFace),
}

struct Scene {
    /// The smallest cell width (the uniform spacing on a uniform grid).
    dx: f64,
    grid: Grid,
    fluid: FluidProperties,
    materials: Vec<SolidMaterial>,
    contacts: Vec<ContactResistance>,
    solids: Vec<(u16, Shape)>,
    sources: Vec<Source>,
    faces: [FaceRule; 6],
    gravity: Option<[f64; 3]>,
    expansion: Option<f64>,
    reference: Option<f64>,
    tolerance: f64,
    turbulence: Turbulence,
    /// Named fans inside the domain.
    internal_fans: Vec<(String, InternalFan)>,
    /// Grilles and porous blocks.
    resistances: Vec<FlowResistance>,
    /// Named two-resistor compact components.
    components: Vec<(String, CompactComponent)>,
    max_iterations: usize,
    wall_seconds: f64,
    transient: Option<Transient>,
    /// Emissivity per material (all zero: no radiation).
    emissivity: Vec<f64>,
    /// Surface-to-surface exchange (else escape to the surroundings only)
    /// and its patch tile edge, in faces.
    surface_exchange: bool,
    patch_size: usize,
    /// Rays per exposed face and stream seed of the radiation estimate.
    rays_per_face: usize,
    ray_seed: u64,
}

/// A backward-Euler march of the energy equation over the converged steady
/// flow, with every source scaled by a piecewise-linear schedule.
struct Transient {
    time_step_s: f64,
    steps: usize,
    initial_temperature_k: Option<f64>,
    /// `(time s, power scale)`, increasing time; constant beyond the ends.
    schedule: Vec<(f64, f64)>,
    /// `flow: "unsteady"`: march the flow with the energy (else the energy
    /// marches over the steady flow).
    unsteady: Option<UnsteadyOptions>,
}

/// Controls of the unsteady flow march.
struct UnsteadyOptions {
    scheme: TimeScheme,
    inner_iterations: usize,
    inner_tolerance: f64,
    /// `(time s, inlet velocity scale)`, like the power schedule.
    inlet_schedule: Vec<(f64, f64)>,
    /// `energy: "steady-on-mean-flow"`: march the flow only, then solve the
    /// steady energy equation on its time-averaged fluxes.
    mean_flow_energy: bool,
}

/// A piecewise-linear `[[time_s, scale], ...]` schedule (increasing times).
fn parse_schedule(value: Option<&J>, key: &str) -> Result<Vec<(f64, f64)>> {
    let mut schedule = Vec::new();
    if let Some(points) = value {
        for point in points
            .as_array()
            .ok_or_else(|| bad(format!("transient.{key} must be [[time_s, scale], ...]")))?
        {
            let pair = point
                .as_array()
                .filter(|pair| pair.len() == 2)
                .and_then(|pair| Some((pair[0].as_f64()?, pair[1].as_f64()?)))
                .filter(|(t, k)| t.is_finite() && k.is_finite())
                .ok_or_else(|| bad(format!("transient.{key} entries must be [time_s, scale]")))?;
            if schedule
                .last()
                .is_some_and(|&(t, _): &(f64, f64)| pair.0 <= t)
            {
                return Err(bad(format!("transient.{key} times must increase")));
            }
            schedule.push(pair);
        }
    }
    Ok(schedule)
}

/// Linear interpolation in `schedule`, constant beyond the ends (1 when
/// empty).
fn interpolate(schedule: &[(f64, f64)], time: f64) -> f64 {
    let Some(&(t0, k0)) = schedule.first() else {
        return 1.0;
    };
    if time <= t0 {
        return k0;
    }
    for pair in schedule.windows(2) {
        let ((ta, ka), (tb, kb)) = (pair[0], pair[1]);
        if time <= tb {
            return ka + (kb - ka) * (time - ta) / (tb - ta);
        }
    }
    schedule.last().map_or(1.0, |&(_, k)| k)
}

impl Transient {
    const MAX_STEPS: usize = 100_000;

    fn parse(value: &J) -> Result<Self> {
        let at = "transient";
        let time_step_s = number(value, "time_step_s", at)?;
        let steps = number(value, "steps", at)?;
        if time_step_s <= 0.0
            || steps < 1.0
            || steps > Self::MAX_STEPS as f64
            || steps.fract() != 0.0
        {
            return Err(bad(format!(
                "transient needs time_step_s > 0 and a whole number of steps in 1..={}",
                Self::MAX_STEPS
            )));
        }
        let schedule = parse_schedule(value.get("power_schedule"), "power_schedule")?;
        let unsteady = match value.str_field("flow") {
            None | Some("frozen") => None,
            Some("unsteady") => {
                let scheme = match value.str_field("scheme") {
                    None | Some("bdf2") => TimeScheme::Bdf2,
                    Some("backward-euler") => TimeScheme::BackwardEuler,
                    Some(other) => {
                        return Err(bad(format!(
                            "transient.scheme must be \"bdf2\" or \"backward-euler\", not {other}"
                        )));
                    }
                };
                let inner = optional_number(value, "inner_iterations", at)?.unwrap_or(100.0);
                let tolerance = optional_number(value, "inner_tolerance", at)?.unwrap_or(1e-6);
                if !(1.0..=100_000.0).contains(&inner) || inner.fract() != 0.0 || !(tolerance > 0.0)
                {
                    return Err(bad(
                        "transient.inner_iterations must be a whole number >= 1 and inner_tolerance positive",
                    ));
                }
                let mean_flow_energy = match value.str_field("energy") {
                    None | Some("march") => false,
                    Some("steady-on-mean-flow") => true,
                    Some(other) => {
                        return Err(bad(format!(
                            "transient.energy must be \"march\" or \"steady-on-mean-flow\", not {other}"
                        )));
                    }
                };
                Some(UnsteadyOptions {
                    scheme,
                    inner_iterations: inner as usize,
                    inner_tolerance: tolerance,
                    inlet_schedule: parse_schedule(value.get("inlet_schedule"), "inlet_schedule")?,
                    mean_flow_energy,
                })
            }
            Some(other) => {
                return Err(bad(format!(
                    "transient.flow must be \"frozen\" or \"unsteady\", not {other}"
                )));
            }
        };
        Ok(Self {
            time_step_s,
            steps: steps as usize,
            initial_temperature_k: optional_number(value, "initial_temperature_k", at)?,
            schedule,
            unsteady,
        })
    }

    fn scale(&self, time: f64) -> f64 {
        interpolate(&self.schedule, time)
    }
}

const FACE_KEYS: [&str; 6] = ["x-", "x+", "y-", "y+", "z-", "z+"];
/// The hidden material of two-resistor component boxes (they block flow;
/// their conduction is the compact model's).
const COMPACT_MATERIAL: &str = "compact-model";

impl Scene {
    #[allow(clippy::too_many_lines)] // one schema, field by field
    fn from_root(root: &J, base: &std::path::Path) -> Result<Self> {
        let root = root.clone();
        if root.str_field("schema") != Some(SCHEMA) {
            return Err(bad(format!("schema must be {SCHEMA}")));
        }
        let grid = match root.get("grid") {
            Some(zones) => {
                if root.get("voxel_m").is_some() {
                    return Err(bad(
                        "declare either voxel_m (uniform) or grid (graded), not both",
                    ));
                }
                parse_grid(zones)?
            }
            None => {
                let dx = number(&root, "voxel_m", "scene")?;
                if dx <= 0.0 {
                    return Err(bad("scene.voxel_m must be positive"));
                }
                let size = vec3(&root, "size_m", "scene")?;
                let mut widths = [Vec::new(), Vec::new(), Vec::new()];
                for a in 0..3 {
                    let cells = size[a] / dx;
                    let rounded = cells.round();
                    if rounded < 1.0 || (cells - rounded).abs() > 1e-6 * rounded.max(1.0) {
                        return Err(bad(format!(
                            "scene.size_m[{a}] = {} is not a whole number of {dx} m voxels",
                            size[a]
                        )));
                    }
                    widths[a] = vec![dx; rounded as usize];
                }
                Grid::from_widths(widths)
            }
        };
        let dims = grid.dims();
        let dx = grid.min_width;
        if dims.iter().product::<usize>() > MAX_CELLS {
            return Err(Failure {
                code: "cooling-cht-budget",
                message: format!(
                    "{} cells exceed the {MAX_CELLS}-cell cap",
                    dims.iter().product::<usize>()
                ),
            });
        }
        let fluid = match root.get("fluid") {
            None => FluidProperties::dry_air_300k(),
            Some(value) if value.as_str() == Some("dry-air-300k") => {
                FluidProperties::dry_air_300k()
            }
            Some(value) if value.as_object().is_some() => FluidProperties {
                density_kg_m3: number(value, "density_kg_m3", "fluid")?,
                specific_heat_j_kg_k: number(value, "specific_heat_j_kg_k", "fluid")?,
                conductivity_w_m_k: number(value, "conductivity_w_m_k", "fluid")?,
                kinematic_viscosity_m2_s: number(value, "kinematic_viscosity_m2_s", "fluid")?,
            },
            Some(_) => {
                return Err(bad(
                    "fluid must be \"dry-air-300k\" or an object of properties",
                ));
            }
        };
        let mut materials = Vec::new();
        let mut emissivities = Vec::new();
        let mut names = Vec::new();
        for (i, item) in root
            .get("materials")
            .and_then(J::as_array)
            .unwrap_or(&[])
            .iter()
            .enumerate()
        {
            let at = format!("materials[{i}]");
            let name = item
                .str_field("name")
                .ok_or_else(|| bad(format!("{at}.name is required")))?;
            if names.contains(&name) {
                return Err(bad(format!("{at}: duplicate material name {name}")));
            }
            if name == COMPACT_MATERIAL {
                return Err(bad(format!("{at}: {COMPACT_MATERIAL} is reserved")));
            }
            names.push(name);
            // A number (isotropic) or [k_x, k_y, k_z] (orthotropic along the
            // grid axes, e.g. a PCB laminate).
            let material = match item.get("conductivity_w_m_k").and_then(J::as_array) {
                Some(_) => {
                    let k = vec3(item, "conductivity_w_m_k", &at)?;
                    SolidMaterial::new(name, (k[0] * k[1] * k[2]).cbrt()).with_orthotropic(k)
                }
                None => SolidMaterial::new(name, number(item, "conductivity_w_m_k", &at)?),
            };
            let material = match optional_number(item, "volumetric_heat_capacity_j_m3_k", &at)? {
                Some(rho_c) => material.with_heat_capacity(rho_c),
                None => material,
            };
            let emissivity = optional_number(item, "emissivity", &at)?.unwrap_or(0.0);
            if !(0.0..=1.0).contains(&emissivity) {
                return Err(bad(format!("{at}.emissivity must lie in [0, 1]")));
            }
            emissivities.push(emissivity);
            materials.push(material);
        }
        if materials.len() > usize::from(u16::MAX) {
            return Err(bad("at most 65535 materials"));
        }
        let mut contacts = Vec::new();
        for (i, item) in root
            .get("contacts")
            .and_then(J::as_array)
            .unwrap_or(&[])
            .iter()
            .enumerate()
        {
            let at = format!("contacts[{i}]");
            let pair = item
                .get("between")
                .and_then(J::as_array)
                .filter(|pair| pair.len() == 2)
                .ok_or_else(|| bad(format!("{at}.between must name two materials")))?;
            let mut index = [0u16; 2];
            for (slot, entry) in index.iter_mut().zip(pair) {
                let name = entry
                    .as_str()
                    .ok_or_else(|| bad(format!("{at}.between entries must be material names")))?;
                *slot = names
                    .iter()
                    .position(|n| *n == name)
                    .and_then(|p| u16::try_from(p).ok())
                    .ok_or_else(|| bad(format!("{at}: unknown material {name}")))?;
            }
            if index[0] == index[1] {
                return Err(bad(format!(
                    "{at}: a contact joins two different materials"
                )));
            }
            let resistance = number(item, "resistance_m2_k_w", &at)?;
            if resistance < 0.0 {
                return Err(bad(format!("{at}.resistance_m2_k_w must be non-negative")));
            }
            contacts.push(ContactResistance {
                materials: (index[0], index[1]),
                resistance_m2_k_w: resistance,
            });
        }
        let mut solids = Vec::new();
        for (i, item) in root
            .get("solids")
            .and_then(J::as_array)
            .unwrap_or(&[])
            .iter()
            .enumerate()
        {
            let at = format!("solids[{i}]");
            let material = item
                .str_field("material")
                .ok_or_else(|| bad(format!("{at}.material is required")))?;
            let index = names
                .iter()
                .position(|n| *n == material)
                .ok_or_else(|| bad(format!("{at}: unknown material {material}")))?;
            let material = u16::try_from(index).expect("bounded above");
            if let Some(heatsink) = item.get("heatsink") {
                for part in plate_fin_heatsink(heatsink, &grid, &format!("{at}.heatsink"))? {
                    solids.push((material, Shape::Box(part)));
                }
                continue;
            }
            let shape = if item.get("stl").is_some() {
                Shape::Mesh(Box::new(MeshSolid::load(item, &at, base)?))
            } else {
                Shape::Box(Aabb::parse(item, &at)?)
            };
            solids.push((material, shape));
        }
        // Two-resistor components: voxelized after (over) the declared
        // solids with the hidden compact material.
        let mut components = Vec::new();
        let component_items = array_of(&root, "components")?;
        if !component_items.is_empty() {
            let index = u16::try_from(materials.len()).map_err(|_| bad("too many materials"))?;
            materials.push(SolidMaterial::new(COMPACT_MATERIAL, 1.0));
            emissivities.push(0.0);
            for (i, item) in component_items.iter().enumerate() {
                let at = format!("components[{i}]");
                let region = Aabb::parse(item, &at)?;
                let (lo, hi) = cell_box(&region, &grid, &at)?;
                let board_face = item
                    .str_field("board_side")
                    .and_then(|key| FACE_KEYS.iter().position(|k| *k == key))
                    .map(|side| Face3::ALL[side])
                    .ok_or_else(|| {
                        bad(format!("{at}.board_side must be x-, x+, y-, y+, z- or z+"))
                    })?;
                let power_w = number(item, "power_w", &at)?;
                let case = number(item, "junction_to_case_k_w", &at)?;
                let board = number(item, "junction_to_board_k_w", &at)?;
                if power_w < 0.0 || !(case > 0.0) || !(board > 0.0) {
                    return Err(bad(format!(
                        "{at}: power_w must be non-negative and both resistances positive"
                    )));
                }
                solids.push((index, Shape::Box(region)));
                components.push((
                    item.str_field("name").unwrap_or("component").to_string(),
                    CompactComponent {
                        lo,
                        hi,
                        board_face,
                        power_w,
                        junction_to_case_k_w: case,
                        junction_to_board_k_w: board,
                    },
                ));
            }
        }
        let mut sources = Vec::new();
        for (i, item) in root
            .get("sources")
            .and_then(J::as_array)
            .unwrap_or(&[])
            .iter()
            .enumerate()
        {
            let at = format!("sources[{i}]");
            let power_w = number(item, "power_w", &at)?;
            if power_w < 0.0 {
                return Err(bad(format!("{at}.power_w must be non-negative")));
            }
            sources.push(Source {
                name: item.str_field("name").unwrap_or("source").to_string(),
                power_w,
                region: Aabb::parse(item, &at)?,
            });
        }
        let mut faces = [FaceRule::Wall(ThermalFace::Adiabatic); 6];
        let face_table = root.get("faces");
        if let Some(table) = face_table
            && let Some(entries) = table.as_object()
        {
            for (key, _) in entries {
                if !FACE_KEYS.contains(&key.as_str()) {
                    return Err(bad(format!("faces.{key}: use x-, x+, y-, y+, z-, z+")));
                }
            }
        }
        for (slot, key) in faces.iter_mut().zip(FACE_KEYS) {
            let Some(rule) = face_table.and_then(|t| t.get(key)) else {
                continue;
            };
            let at = format!("faces.{key}");
            *slot = match rule.str_field("type") {
                Some("inlet") => FaceRule::Inlet {
                    velocity: vec3(rule, "velocity_m_s", &at)?,
                    temperature: number(rule, "temperature_k", &at)?,
                },
                Some("fan") => FaceRule::Fan {
                    curve: fan_curve(rule, &at)?,
                    temperature: number(rule, "temperature_k", &at)?,
                },
                Some("opening") => FaceRule::Opening {
                    ambient: number(rule, "ambient_k", &at)?,
                },
                Some("symmetry") => FaceRule::Symmetry,
                Some("wall") => {
                    let temperature = optional_number(rule, "temperature_k", &at)?;
                    let flux = optional_number(rule, "heat_flux_w_m2", &at)?;
                    let htc = optional_number(rule, "htc_w_m2_k", &at)?;
                    FaceRule::Wall(match (temperature, flux, htc) {
                        (None, None, None) => ThermalFace::Adiabatic,
                        (Some(t), None, None) => ThermalFace::Temperature(t),
                        (None, Some(q), None) => ThermalFace::HeatFlux(q),
                        (None, None, Some(h)) => ThermalFace::Convective {
                            h,
                            ambient: number(rule, "ambient_k", &at)?,
                        },
                        _ => {
                            return Err(bad(format!(
                                "{at}: declare at most one of temperature_k, heat_flux_w_m2, htc_w_m2_k"
                            )));
                        }
                    })
                }
                _ => {
                    return Err(bad(format!(
                        "{at}.type must be inlet, fan, opening, symmetry, or wall"
                    )));
                }
            };
        }
        let gravity = match root.get("gravity_m_s2") {
            None => None,
            Some(_) => Some(vec3(&root, "gravity_m_s2", "scene")?),
        };
        let solver = root.get("solver");
        let tolerance = match solver {
            Some(s) => optional_number(s, "tolerance", "solver")?.unwrap_or(1e-6),
            None => 1e-6,
        };
        let max_iterations = match solver {
            Some(s) => optional_number(s, "max_iterations", "solver")?.unwrap_or(5000.0),
            None => 5000.0,
        };
        if !(tolerance > 0.0) || !(max_iterations >= 1.0) {
            return Err(bad(
                "solver.tolerance must be positive and max_iterations at least 1",
            ));
        }
        let turbulence = match solver.and_then(|s| s.str_field("turbulence")) {
            None | Some("laminar") => Turbulence::Laminar,
            Some("lvel") => Turbulence::Lvel,
            Some(other) => {
                return Err(bad(format!(
                    "solver.turbulence must be \"laminar\" or \"lvel\", not {other}"
                )));
            }
        };
        let mut internal_fans = Vec::new();
        for (i, item) in array_of(&root, "internal_fans")?.iter().enumerate() {
            let at = format!("internal_fans[{i}]");
            let patch = face_patch(item, &grid, &at)?;
            let blows_positive = match item.str_field("direction") {
                Some("+") => true,
                Some("-") => false,
                _ => return Err(bad(format!("{at}.direction must be \"+\" or \"-\""))),
            };
            internal_fans.push((
                item.str_field("name").unwrap_or("fan").to_string(),
                InternalFan {
                    patch,
                    blows_positive,
                    curve: fan_curve(item, &at)?,
                },
            ));
        }
        let mut resistances = Vec::new();
        for (i, item) in array_of(&root, "resistances")?.iter().enumerate() {
            let at = format!("resistances[{i}]");
            resistances.push(match item.str_field("type") {
                Some("grille") => {
                    let loss = match (
                        optional_number(item, "loss_coefficient", &at)?,
                        optional_number(item, "free_area_ratio", &at)?,
                    ) {
                        (Some(k), None) if k >= 0.0 => k,
                        (None, Some(f)) => FlowResistance::perforated_plate_loss(f).ok_or_else(
                            || bad(format!("{at}.free_area_ratio must lie in (0, 1]")),
                        )?,
                        _ => {
                            return Err(bad(format!(
                                "{at}: declare one of loss_coefficient (non-negative) or free_area_ratio"
                            )));
                        }
                    };
                    FlowResistance::Planar {
                        patch: face_patch(item, &grid, &at)?,
                        loss_coefficient: loss,
                    }
                }
                Some("porous") => {
                    let (lo, hi) = cell_box(&Aabb::parse(item, &at)?, &grid, &at)?;
                    let permeability_m2 = per_axis(item, "permeability_m2", &at, f64::INFINITY)?;
                    let inertial_per_m = per_axis(item, "inertial_per_m", &at, 0.0)?;
                    if permeability_m2.iter().any(|k| !(*k > 0.0))
                        || inertial_per_m.iter().any(|c| !(*c >= 0.0))
                    {
                        return Err(bad(format!(
                            "{at}: permeability_m2 must be positive and inertial_per_m non-negative"
                        )));
                    }
                    FlowResistance::Volume {
                        lo,
                        hi,
                        permeability_m2,
                        inertial_per_m,
                    }
                }
                _ => return Err(bad(format!("{at}.type must be grille or porous"))),
            });
        }
        let wall_seconds = match root.get("limits") {
            Some(l) => optional_number(l, "wall_seconds", "limits")?.unwrap_or(3600.0),
            None => 3600.0,
        };
        if !(wall_seconds > 0.0) {
            return Err(bad("limits.wall_seconds must be positive"));
        }
        Ok(Self {
            dx,
            grid,
            fluid,
            materials,
            contacts,
            solids,
            sources,
            faces,
            gravity,
            expansion: optional_number(&root, "expansion_per_k", "scene")?,
            reference: optional_number(&root, "reference_temperature_k", "scene")?,
            tolerance,
            turbulence,
            internal_fans,
            resistances,
            components,
            max_iterations: max_iterations as usize,
            wall_seconds,
            transient: root.get("transient").map(Transient::parse).transpose()?,
            emissivity: emissivities,
            surface_exchange: match root
                .get("radiation")
                .and_then(|r| r.get("surface_exchange"))
            {
                None => true,
                Some(value) => match value {
                    J::Bool(flag) => *flag,
                    _ => return Err(bad("radiation.surface_exchange must be true or false")),
                },
            },
            patch_size: match root.get("radiation") {
                Some(r) => {
                    let size = optional_number(r, "patch_size", "radiation")?.unwrap_or(4.0);
                    if !(1.0..=1024.0).contains(&size) || size.fract() != 0.0 {
                        return Err(bad(
                            "radiation.patch_size must be a whole number in 1..=1024",
                        ));
                    }
                    size as usize
                }
                None => 4,
            },
            rays_per_face: match root.get("radiation") {
                Some(r) => {
                    let rays = optional_number(r, "rays_per_face", "radiation")?.unwrap_or(256.0);
                    if !(1.0..=65536.0).contains(&rays) || rays.fract() != 0.0 {
                        return Err(bad(
                            "radiation.rays_per_face must be a whole number in 1..=65536",
                        ));
                    }
                    rays as usize
                }
                None => 256,
            },
            ray_seed: match root.get("radiation") {
                Some(r) => {
                    optional_number(r, "seed", "radiation")?.map_or(0x5EED_0FA1, |v| v as u64)
                }
                None => 0x5EED_0FA1,
            },
        })
    }
}

/// The unsteady march (`transient.flow: "unsteady"`): flow and energy
/// advance together, so there is no steady solution to report; the result
/// carries the march's closure, peak and final temperatures instead.
#[allow(clippy::too_many_lines)] // solve and one linear report
fn execute_unsteady(
    scene: &Scene,
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    flow_config: &SimpleConfig,
    (transient, options): (&Transient, &UnsteadyOptions),
    gate: &CancelGate,
    json_mode: bool,
) -> Result<String> {
    let ambient = scene.faces.iter().find_map(|rule| match rule {
        FaceRule::Inlet { temperature, .. } | FaceRule::Fan { temperature, .. } => {
            Some(*temperature)
        }
        FaceRule::Opening { ambient } => Some(*ambient),
        _ => None,
    });
    let initial = transient
        .initial_temperature_k
        .or(ambient)
        .ok_or_else(|| bad("transient.initial_temperature_k is required here"))?;
    let buoyancy = match scene.gravity.filter(|g| g.iter().any(|v| *v != 0.0)) {
        Some(gravity) => {
            let reference = scene.reference.or(ambient).ok_or_else(|| {
                bad("a buoyant scene needs reference_temperature_k, an inlet, or an opening")
            })?;
            Some(Boussinesq {
                gravity_m_s2: gravity,
                expansion_per_k: scene.expansion.unwrap_or(1.0 / reference),
                reference_temperature_k: reference,
            })
        }
        None => None,
    };
    let mut unsteady = UnsteadyConfig::new(transient.time_step_s, transient.steps);
    unsteady.scheme = options.scheme;
    unsteady.inner_iterations = options.inner_iterations;
    unsteady.inner_tolerance = options.inner_tolerance;
    let started = Instant::now();
    let run = march_conjugate(
        domain,
        &scene.fluid,
        &scene.materials,
        setup,
        flow_config,
        &unsteady,
        buoyancy.as_ref(),
        &vec![initial; domain.cell_count()],
        |t| interpolate(&options.inlet_schedule, t),
        |t| transient.scale(t),
        &EnergyConfig::default(),
        gate,
    )
    .map_err(|e| solver_failure(&e))?;
    let wall_s = started.elapsed().as_secs_f64();
    let records = &run.energy_records;
    let steps = &run.flow.records;
    let worst_closure = records.iter().fold(0.0f64, |m, r| m.max(r.closure_j.abs()));
    let peak = records
        .iter()
        .map(|r| r.max_solid_temperature_k)
        .filter(|t| t.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    let sweeps: usize = steps.iter().map(|r| r.inner_iterations).sum();
    let max_step_sweeps = steps.iter().map(|r| r.inner_iterations).max().unwrap_or(0);
    let last = steps.last().copied();
    let report = &run.flow.flow.report;
    let tie = 1e-6 * scene.dx;
    let solids_max = |field: &[f64]| {
        (0..domain.cell_count())
            .filter(|&c| !domain.is_fluid(c))
            .map(|c| field[c])
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let source_rows: Vec<(String, f64, f64, f64)> = scene
        .sources
        .iter()
        .map(|source| {
            let cells: Vec<usize> = (0..domain.cell_count())
                .filter(|&c| {
                    let [x, y, z] = domain.coords(c);
                    !domain.is_fluid(c)
                        && source.region.contains(probe(domain.center(x, y, z), tie))
                })
                .collect();
            let max = |field: &[f64]| {
                cells
                    .iter()
                    .map(|&c| field[c])
                    .fold(f64::NEG_INFINITY, f64::max)
            };
            (
                source.name.clone(),
                source.power_w,
                max(&run.temperature),
                max(&run.mean_temperature),
            )
        })
        .collect();
    let solver_name = if buoyancy.is_some() {
        "fv-simplec-unsteady-boussinesq"
    } else {
        "fv-simplec-unsteady"
    };
    let scheme = match options.scheme {
        TimeScheme::Bdf2 => "bdf2",
        TimeScheme::BackwardEuler => "backward-euler",
    };
    if json_mode {
        let mut out = format!(
            "{{\"schema\":{},\"status\":\"completed\",\"solver\":{},\"cells\":{},\"fluid_cells\":{},\"voxel_m\":{},\"graded\":{}",
            quote(RESULT_SCHEMA),
            quote(solver_name),
            domain.cell_count(),
            domain.fluid_count(),
            num(scene.dx)?,
            scene.grid.uniform.is_none()
        );
        let _ = write!(
            out,
            ",\"flow\":{{\"steps\":{},\"sweeps\":{sweeps},\"max_step_sweeps\":{max_step_sweeps},\"final_mass_residual\":{},\"final_momentum_residual\":{},\"inflow_m3_s\":{},\"outflow_m3_s\":{},\"final_kinetic_energy_j\":{}}}",
            steps.len(),
            num(last.map_or(0.0, |r| r.mass_residual))?,
            num(last.map_or(0.0, |r| r.momentum_residual))?,
            num(report.inflow_m3_s)?,
            num(report.outflow_m3_s)?,
            num(last.map_or(0.0, |r| r.kinetic_energy_j))?
        );
        let every = records.len().div_ceil(200).max(1);
        let _ = write!(
            out,
            ",\"transient\":{{\"time_step_s\":{},\"steps\":{},\"scheme\":{},\"averaged_steps\":{},\"worst_step_closure_j\":{},\"peak_solid_temperature_k\":{},\"final_max_solid_temperature_k\":{},\"records\":[",
            num(transient.time_step_s)?,
            records.len(),
            quote(scheme),
            run.flow.averaged_steps,
            num(worst_closure)?,
            num_or_null(peak),
            num_or_null(solids_max(&run.temperature))
        );
        let mut first = true;
        for (i, (record, step)) in records.iter().zip(steps).enumerate() {
            if (i + 1) % every != 0 && i + 1 != records.len() {
                continue;
            }
            if !first {
                out.push(',');
            }
            first = false;
            let _ = write!(
                out,
                "{{\"time_s\":{},\"max_solid_temperature_k\":{},\"kinetic_energy_j\":{},\"sweeps\":{}}}",
                num(record.time_s)?,
                num_or_null(record.max_solid_temperature_k),
                num(step.kinetic_energy_j)?,
                step.inner_iterations
            );
        }
        out.push_str("]},\"sources\":[");
        for (i, (name, power, last_max, mean_max)) in source_rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"power_w\":{},\"final_max_temperature_k\":{},\"mean_max_temperature_k\":{}}}",
                quote(name),
                num(*power)?,
                num(*last_max)?,
                num(*mean_max)?
            );
        }
        let _ = writeln!(
            out,
            "],\"max_solid_temperature_k\":{},\"wall_s\":{},\"evidence\":\"Estimated\",\"no_claim\":{}}}",
            num_or_null(solids_max(&run.temperature)),
            num(wall_s)?,
            quote(NO_CLAIM)
        );
        Ok(out)
    } else {
        let mut out = format!(
            "status=completed\nsolver={solver_name}\ncells={}\nfluid_cells={}\nsteps={}\nsweeps={sweeps}\nscheme={scheme}\nworst_step_closure_j={worst_closure:e}\npeak_solid_temperature_k={peak}\nfinal_max_solid_temperature_k={}\n",
            domain.cell_count(),
            domain.fluid_count(),
            records.len(),
            solids_max(&run.temperature)
        );
        for (name, power, last_max, mean_max) in &source_rows {
            let _ = writeln!(
                out,
                "source={name} power_w={power} final_max_temperature_k={last_max:.4} mean_max_temperature_k={mean_max:.4}"
            );
        }
        let _ = writeln!(
            out,
            "wall_s={wall_s:.3}\nevidence=Estimated\nno_claim={NO_CLAIM}"
        );
        Ok(out)
    }
}

#[allow(clippy::too_many_lines)] // build, solve, and one linear report
fn execute(scene: &Scene, gate: &CancelGate, json_mode: bool) -> Result<String> {
    let tie = 1e-6 * scene.dx;
    let domain = scene
        .grid
        .domain(|p| {
            scene
                .solids
                .iter()
                .rev()
                .find(|(_, shape)| shape.contains(probe(p, tie)))
                .map_or(Voxel::Fluid, |(material, _)| Voxel::Solid(*material))
        })
        .map_err(|e| solver_failure(&e))?;
    let fans: Vec<usize> = (0..6)
        .filter(|&side| matches!(scene.faces[side], FaceRule::Fan { .. }))
        .collect();
    if fans.len() > 1 {
        return Err(bad("at most one fan face"));
    }
    let flow_faces: [FvBoundary; 6] = std::array::from_fn(|side| match scene.faces[side] {
        FaceRule::Inlet { velocity, .. } => FvBoundary::Inlet { velocity },
        FaceRule::Fan { curve, .. } => {
            // Initial guess: half the free delivery through the whole face.
            let axis = side / 2;
            let area: f64 = (0..3)
                .filter(|&a| a != axis)
                .map(|a| scene.grid.length(a))
                .product();
            let free_delivery = -curve.pressure(0.0) / curve.slope(0.0);
            let speed = 0.5 * free_delivery.max(0.0) / area;
            let mut velocity = [0.0; 3];
            velocity[axis] = if side % 2 == 0 { speed } else { -speed };
            FvBoundary::Inlet { velocity }
        }
        FaceRule::Opening { .. } => FvBoundary::Outlet,
        FaceRule::Symmetry => FvBoundary::Symmetry,
        FaceRule::Wall(_) => FvBoundary::wall(),
    });
    let thermal_faces = scene.faces.map(|rule| match rule {
        FaceRule::Inlet { temperature, .. } | FaceRule::Fan { temperature, .. } => {
            ThermalFace::Inflow { temperature }
        }
        FaceRule::Opening { ambient } => ThermalFace::Outflow {
            backflow_temperature: ambient,
        },
        FaceRule::Symmetry => ThermalFace::Adiabatic,
        FaceRule::Wall(thermal) => thermal,
    });
    let mut setup = ThermalSetup::new(thermal_faces);
    setup.contacts.clone_from(&scene.contacts);
    setup.compact_components = scene.components.iter().map(|(_, part)| *part).collect();
    let mut source_cells = Vec::with_capacity(scene.sources.len());
    for source in &scene.sources {
        let count = setup.add_uniform_power(&domain, source.power_w, |p| {
            source.region.contains(probe(p, tie))
                && !domain.is_fluid({
                    let cell = |a: usize| scene.grid.locate(a, p[a]);
                    domain.index(cell(0), cell(1), cell(2))
                })
        });
        if count == 0 {
            return Err(bad(format!(
                "source {} covers no solid voxel centre",
                source.name
            )));
        }
        source_cells.push(count);
    }
    let mut flow_config = SimpleConfig::new(flow_faces);
    flow_config.tolerance = scene.tolerance;
    flow_config.turbulence = scene.turbulence;
    flow_config.max_iterations = scene.max_iterations;
    flow_config.internal_fans = scene.internal_fans.iter().map(|(_, fan)| *fan).collect();
    flow_config.resistances.clone_from(&scene.resistances);
    flow_config.fan = fans.first().map(|&side| FanInlet {
        face: Face3::ALL[side],
        curve: match scene.faces[side] {
            FaceRule::Fan { curve, .. } => curve,
            _ => unreachable!("filtered to fan faces"),
        },
    });
    // Surface radiation: openings, inlets and fans are the surroundings at
    // their temperatures; walls and symmetry planes are opaque.
    let radiation = scene.emissivity.iter().any(|e| *e > 0.0).then(|| {
        let surroundings = scene.faces.map(|rule| match rule {
            FaceRule::Inlet { temperature, .. } | FaceRule::Fan { temperature, .. } => {
                Some(temperature)
            }
            FaceRule::Opening { ambient } => Some(ambient),
            FaceRule::Symmetry | FaceRule::Wall(_) => None,
        });
        let mut config = RadiationConfig::new(scene.emissivity.clone(), surroundings);
        config.rays_per_face = scene.rays_per_face;
        config.seed = scene.ray_seed;
        config.surface_exchange = scene.surface_exchange;
        config.patch_size = scene.patch_size;
        config
    });
    // Steady energy on the time-averaged flow of an unsteady march.
    let mean_flow = scene
        .transient
        .as_ref()
        .and_then(|t| t.unsteady.as_ref().map(|o| (t, o)))
        .filter(|(_, o)| o.mean_flow_energy);
    if radiation.is_some() && scene.transient.is_some() && mean_flow.is_none() {
        return Err(bad(
            "transient marches do not carry radiation; drop the transient block or the emissivities",
        ));
    }
    if let Some(transient) = &scene.transient
        && let Some(options) = &transient.unsteady
        && !options.mean_flow_energy
    {
        return execute_unsteady(
            scene,
            &domain,
            &setup,
            &flow_config,
            (transient, options),
            gate,
            json_mode,
        );
    }
    let mut radiated: Option<(f64, usize, usize)> = None;
    let mut averaged_steps: Option<usize> = None;
    if mean_flow.is_some() && scene.gravity.is_some_and(|g| g.iter().any(|v| *v != 0.0)) {
        return Err(bad(
            "steady-on-mean-flow energy needs a forced flow: a buoyant flow depends on the temperature it would freeze",
        ));
    }
    let started = Instant::now();
    let buoyant = scene.gravity.filter(|g| g.iter().any(|v| *v != 0.0));
    let (flow, energy, couplings): (FvFlow, EnergySolution, Option<usize>) =
        if let Some(gravity) = buoyant {
            let reference = scene
                .reference
                .or_else(|| {
                    scene.faces.iter().find_map(|rule| match rule {
                        FaceRule::Inlet { temperature, .. } | FaceRule::Fan { temperature, .. } => {
                            Some(*temperature)
                        }
                        FaceRule::Opening { ambient } => Some(*ambient),
                        _ => None,
                    })
                })
                .ok_or_else(|| {
                    bad("a buoyant scene needs reference_temperature_k, an inlet, or an opening")
                })?;
            let expansion = scene.expansion.unwrap_or(1.0 / reference);
            let mut config = FvBuoyancyConfig::new(gravity, expansion, reference, flow_config);
            config.max_couplings = scene.max_iterations;
            config.radiation.clone_from(&radiation);
            let run = fv_natural_convection(
                &domain,
                &scene.fluid,
                &scene.materials,
                &setup,
                &config,
                gate,
            )
            .map_err(|e| solver_failure(&e))?;
            if let Some(w) = run.report.radiated_w {
                radiated = Some((w, run.report.couplings, 0));
            }
            (run.flow, run.energy, Some(run.report.couplings))
        } else {
            let flow = match mean_flow {
                Some((transient, options)) => {
                    let mut unsteady = UnsteadyConfig::new(transient.time_step_s, transient.steps);
                    unsteady.scheme = options.scheme;
                    unsteady.inner_iterations = options.inner_iterations;
                    unsteady.inner_tolerance = options.inner_tolerance;
                    let run = simple_unsteady(
                        &domain,
                        &scene.fluid,
                        &flow_config,
                        &unsteady,
                        |t| interpolate(&options.inlet_schedule, t),
                        gate,
                    )
                    .map_err(|e| solver_failure(&e))?;
                    averaged_steps = Some(run.averaged_steps);
                    FvFlow {
                        field: run.mean_field,
                        velocity_m_s: run.mean_velocity_m_s,
                        ..run.flow
                    }
                }
                None => simple_flow(&domain, &scene.fluid, &flow_config, gate)
                    .map_err(|e| solver_failure(&e))?,
            };
            if scene.turbulence != Turbulence::Laminar {
                setup.eddy_conductivity_w_m_k = flow.eddy_conductivity(&scene.fluid);
            }
            let energy = match &radiation {
                Some(config) => {
                    let (energy, report) = solve_energy_radiating(
                        &domain,
                        &scene.fluid,
                        &scene.materials,
                        &flow.field,
                        &setup,
                        &EnergyConfig::default(),
                        config,
                        gate,
                    )
                    .map_err(|e| solver_failure(&e))?;
                    radiated = Some((report.radiated_w, report.iterations, report.exposed_faces));
                    energy
                }
                None => solve_energy(
                    &domain,
                    &scene.fluid,
                    &scene.materials,
                    &flow.field,
                    &setup,
                    &EnergyConfig::default(),
                    gate,
                )
                .map_err(|e| solver_failure(&e))?,
            };
            (flow, energy, None)
        };
    // Optional transient march over the converged (forced) flow.
    let march = match (&scene.transient, couplings) {
        (None, _) => None,
        (Some(_), _) if mean_flow.is_some() => None,
        (Some(_), Some(_)) => {
            return Err(bad(
                "transient runs freeze the flow; a buoyant scene's flow depends on temperature",
            ));
        }
        (Some(transient), None) => {
            let initial = transient
                .initial_temperature_k
                .or_else(|| {
                    scene.faces.iter().find_map(|rule| match rule {
                        FaceRule::Inlet { temperature, .. } | FaceRule::Fan { temperature, .. } => {
                            Some(*temperature)
                        }
                        FaceRule::Opening { ambient } => Some(*ambient),
                        _ => None,
                    })
                })
                .ok_or_else(|| bad("transient.initial_temperature_k is required here"))?;
            let config = TransientConfig {
                time_step_s: transient.time_step_s,
                steps: transient.steps,
                energy: EnergyConfig::default(),
            };
            let solution = march_energy(
                &domain,
                &scene.fluid,
                &scene.materials,
                &flow.field,
                &setup,
                &vec![initial; domain.cell_count()],
                |t| transient.scale(t),
                &config,
                gate,
            )
            .map_err(|e| solver_failure(&e))?;
            Some(solution)
        }
    };
    let wall_s = started.elapsed().as_secs_f64();
    let temperature = &energy.temperature;
    // Per-material and per-source temperatures.
    let mut material_rows = Vec::new();
    for (index, material) in scene.materials.iter().enumerate() {
        if material.label == COMPACT_MATERIAL {
            continue;
        }
        let cells: Vec<usize> = (0..domain.cell_count())
            .filter(|&c| domain.voxel_at(c) == Voxel::Solid(u16::try_from(index).expect("bounded")))
            .collect();
        if cells.is_empty() {
            continue;
        }
        let max = cells
            .iter()
            .map(|&c| temperature[c])
            .fold(f64::NEG_INFINITY, f64::max);
        let mean = cells.iter().map(|&c| temperature[c]).sum::<f64>() / cells.len() as f64;
        material_rows.push((material.label.clone(), cells.len(), max, mean));
    }
    let mut source_rows = Vec::new();
    for (source, count) in scene.sources.iter().zip(&source_cells) {
        let max = (0..domain.cell_count())
            .filter(|&c| {
                let [x, y, z] = domain.coords(c);
                !domain.is_fluid(c) && source.region.contains(probe(domain.center(x, y, z), tie))
            })
            .map(|c| temperature[c])
            .fold(f64::NEG_INFINITY, f64::max);
        source_rows.push((source.name.clone(), source.power_w, *count, max));
    }
    let hottest = (0..domain.cell_count())
        .filter(|&c| !domain.is_fluid(c))
        .max_by(|&a, &b| temperature[a].total_cmp(&temperature[b]))
        .map(|c| {
            let [x, y, z] = domain.coords(c);
            (temperature[c], domain.center(x, y, z))
        });
    let max_speed = flow
        .velocity_m_s
        .iter()
        .map(|v| fs_math::det::sqrt(v.iter().map(|x| x * x).sum::<f64>()))
        .fold(0.0f64, f64::max);
    let r = &flow.report;
    let b = &energy.report.balance;
    let nu = scene.fluid.kinematic_viscosity_m2_s;
    let max_eddy_ratio = flow
        .eddy_viscosity_m2_s
        .iter()
        .fold(0.0f64, |m, nu_t| m.max(nu_t / nu));
    let turbulence_name = match scene.turbulence {
        Turbulence::Laminar => "laminar",
        Turbulence::Lvel => "lvel",
    };
    let solver_name = if couplings.is_some() {
        "fv-simplec-boussinesq"
    } else if averaged_steps.is_some() {
        "fv-simplec-unsteady-mean-flow"
    } else {
        "fv-simplec"
    };
    if json_mode {
        let mut out = format!(
            "{{\"schema\":{},\"status\":\"completed\",\"solver\":{},\"cells\":{},\"fluid_cells\":{},\"voxel_m\":{},\"graded\":{}",
            quote(RESULT_SCHEMA),
            quote(solver_name),
            domain.cell_count(),
            domain.fluid_count(),
            num(scene.dx)?,
            scene.grid.uniform.is_none()
        );
        let _ = write!(
            out,
            ",\"flow\":{{\"iterations\":{},\"mass_residual\":{},\"momentum_residual\":{},\"inflow_m3_s\":{},\"outflow_m3_s\":{},\"max_divergence_m3_s\":{},\"max_speed_m_s\":{},\"max_cell_reynolds\":{}",
            r.iterations,
            num(r.mass_residual)?,
            num(r.momentum_residual)?,
            num(r.inflow_m3_s)?,
            num(r.outflow_m3_s)?,
            num(r.max_divergence_m3_s)?,
            num(max_speed)?,
            num(r.max_cell_reynolds)?
        );
        let _ = write!(out, ",\"turbulence\":{}", quote(turbulence_name));
        if scene.turbulence != Turbulence::Laminar {
            let _ = write!(
                out,
                ",\"max_eddy_viscosity_ratio\":{}",
                num(max_eddy_ratio)?
            );
        }
        if let Some(couplings) = couplings {
            let _ = write!(out, ",\"energy_couplings\":{couplings}");
        }
        if let Some(steps) = averaged_steps {
            let _ = write!(out, ",\"averaged_steps\":{steps}");
        }
        if let Some((q, dp, mismatch)) = r.fan {
            let _ = write!(
                out,
                ",\"fan_flow_m3_s\":{},\"fan_pressure_pa\":{},\"fan_residual\":{}",
                num(q)?,
                num(dp)?,
                num(mismatch)?
            );
        }
        if !r.internal_fans.is_empty() {
            out.push_str(",\"internal_fans\":[");
            for (i, ((name, _), &(q, dp))) in
                scene.internal_fans.iter().zip(&r.internal_fans).enumerate()
            {
                if i > 0 {
                    out.push(',');
                }
                let _ = write!(
                    out,
                    "{{\"name\":{},\"flow_m3_s\":{},\"pressure_rise_pa\":{}}}",
                    quote(name),
                    num(q)?,
                    num(dp)?
                );
            }
            out.push(']');
        }
        let _ = write!(
            out,
            "}},\"energy\":{{\"iterations\":{},\"relative_residual\":{},\"source_w\":{},\"boundary_outflow_w\":{},\"advective_outflow_w\":{},\"sink_outflow_w\":{},\"balance_relative_residual\":{},\"max_cell_peclet\":{}}}",
            energy.report.iterations,
            num(energy.report.relative_residual)?,
            num(b.source_w)?,
            num(b.boundary_outflow_w)?,
            num(b.advective_outflow_w)?,
            num(b.sink_outflow_w)?,
            num(b.relative_residual)?,
            num(energy.report.max_cell_peclet)?
        );
        if let Some((w, iterations, faces)) = radiated {
            let _ = write!(
                out,
                ",\"radiation\":{{\"radiated_w\":{},\"iterations\":{iterations},\"exposed_faces\":{faces},\"rays_per_face\":{},\"surface_exchange\":{}}}",
                num(w)?,
                scene.rays_per_face,
                scene.surface_exchange
            );
        }
        out.push_str(",\"materials\":[");
        for (i, (name, cells, max, mean)) in material_rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"cells\":{cells},\"max_temperature_k\":{},\"mean_temperature_k\":{}}}",
                quote(name),
                num(*max)?,
                num(*mean)?
            );
        }
        out.push_str("],\"sources\":[");
        for (i, (name, power, cells, max)) in source_rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"power_w\":{},\"cells\":{cells},\"max_temperature_k\":{}}}",
                quote(name),
                num(*power)?,
                num(*max)?
            );
        }
        out.push(']');
        if !scene.components.is_empty() {
            out.push_str(",\"components\":[");
            for (i, ((name, part), junction)) in
                scene.components.iter().zip(&energy.junctions).enumerate()
            {
                if i > 0 {
                    out.push(',');
                }
                let _ = write!(
                    out,
                    "{{\"name\":{},\"power_w\":{},\"junction_temperature_k\":{},\"case_w\":{},\"board_w\":{}}}",
                    quote(name),
                    num(part.power_w)?,
                    num(junction.temperature_k)?,
                    num(junction.case_w)?,
                    num(junction.board_w)?
                );
            }
            out.push(']');
        }
        if let Some(march) = &march {
            // At most ~200 records: every k-th step and the last.
            let every = march.records.len().div_ceil(200).max(1);
            let worst_closure = march
                .records
                .iter()
                .fold(0.0f64, |m, r| m.max(r.closure_j.abs()));
            let _ = write!(
                out,
                ",\"transient\":{{\"steps\":{},\"worst_step_closure_j\":{},\"final_max_solid_temperature_k\":{},\"records\":[",
                march.records.len(),
                num(worst_closure)?,
                num(march
                    .records
                    .last()
                    .map_or(f64::NAN, |r| r.max_solid_temperature_k))?
            );
            let mut first = true;
            for (i, record) in march.records.iter().enumerate() {
                if (i + 1) % every != 0 && i + 1 != march.records.len() {
                    continue;
                }
                if !first {
                    out.push(',');
                }
                first = false;
                let _ = write!(
                    out,
                    "{{\"time_s\":{},\"max_solid_temperature_k\":{},\"source_j\":{},\"stored_j\":{}}}",
                    num(record.time_s)?,
                    num(record.max_solid_temperature_k)?,
                    num(record.source_j)?,
                    num(record.stored_energy_change_j)?
                );
            }
            out.push_str("]}");
        }
        if let Some((t, at)) = hottest {
            let _ = write!(
                out,
                ",\"max_solid_temperature_k\":{},\"max_solid_at_m\":[{},{},{}]",
                num(t)?,
                num(at[0])?,
                num(at[1])?,
                num(at[2])?
            );
        }
        let _ = writeln!(
            out,
            ",\"wall_s\":{},\"evidence\":\"Estimated\",\"no_claim\":{}}}",
            num(wall_s)?,
            quote(NO_CLAIM)
        );
        Ok(out)
    } else {
        let mut out = format!(
            "status=completed\nsolver={solver_name}\ncells={}\nfluid_cells={}\nflow_iterations={}\nmass_residual={:e}\nmomentum_residual={:e}\ninflow_m3_s={:e}\noutflow_m3_s={:e}\nmax_speed_m_s={max_speed}\nenergy_balance_relative_residual={:e}\n",
            domain.cell_count(),
            domain.fluid_count(),
            r.iterations,
            r.mass_residual,
            r.momentum_residual,
            r.inflow_m3_s,
            r.outflow_m3_s,
            b.relative_residual
        );
        if scene.turbulence != Turbulence::Laminar {
            let _ = writeln!(
                out,
                "turbulence={turbulence_name}\nmax_eddy_viscosity_ratio={max_eddy_ratio}"
            );
        }
        if let Some((q, dp, _)) = r.fan {
            let _ = writeln!(out, "fan_flow_m3_s={q:e}\nfan_pressure_pa={dp}");
        }
        for ((name, _), (q, dp)) in scene.internal_fans.iter().zip(&r.internal_fans) {
            let _ = writeln!(
                out,
                "internal_fan.{name}.flow_m3_s={q:e}\ninternal_fan.{name}.pressure_rise_pa={dp}"
            );
        }
        if let Some(record) = march.as_ref().and_then(|m| m.records.last()) {
            let _ = writeln!(
                out,
                "transient_steps={}\ntransient_final_time_s={}\ntransient_final_max_solid_temperature_k={:.4}",
                march.as_ref().map_or(0, |m| m.records.len()),
                record.time_s,
                record.max_solid_temperature_k
            );
        }
        for (name, cells, max, mean) in &material_rows {
            let _ = writeln!(
                out,
                "material={name} cells={cells} max_temperature_k={max:.4} mean_temperature_k={mean:.4}"
            );
        }
        for (name, power, cells, max) in &source_rows {
            let _ = writeln!(
                out,
                "source={name} power_w={power} cells={cells} max_temperature_k={max:.4}"
            );
        }
        for ((name, part), junction) in scene.components.iter().zip(&energy.junctions) {
            let _ = writeln!(
                out,
                "component={name} power_w={} junction_temperature_k={:.4} case_w={:.6} board_w={:.6}",
                part.power_w, junction.temperature_k, junction.case_w, junction.board_w
            );
        }
        if let Some((t, at)) = hottest {
            let _ = writeln!(
                out,
                "max_solid_temperature_k={t:.4} at_m=({:.6},{:.6},{:.6})",
                at[0], at[1], at[2]
            );
        }
        let _ = writeln!(
            out,
            "wall_s={wall_s:.2}\nevidence=Estimated\nno_claim={NO_CLAIM}"
        );
        Ok(out)
    }
}

fn diagnostic(code: u8, failure: &Failure, json_mode: bool) -> CommandOutput {
    let stderr = if json_mode {
        format!(
            "{{\"schema\":\"frankensim.cooling-cht.diagnostic.v1\",\"code\":{},\"message\":{}}}\n",
            quote(failure.code),
            quote(&failure.message)
        )
    } else {
        format!("{failure}\n")
    };
    CommandOutput {
        exit_code: code,
        stdout: String::new(),
        stderr,
    }
}

/// Execute `scene` under its wall-time limit: the limit trips the same
/// cancellation gate the solvers poll, so an exhausted budget publishes
/// nothing.
fn execute_within_budget(scene: &Scene, json_mode: bool) -> Result<String> {
    let gate = CancelGate::new();
    let duration = Duration::from_secs_f64(scene.wall_seconds);
    std::thread::scope(|scope| {
        let (stop, stopped) = std::sync::mpsc::channel::<()>();
        let gate_ref = &gate;
        scope.spawn(move || {
            if matches!(
                stopped.recv_timeout(duration),
                Err(std::sync::mpsc::RecvTimeoutError::Timeout)
            ) {
                gate_ref.request();
            }
        });
        let result = execute(scene, &gate, json_mode);
        let _ = stop.send(());
        result
    })
}

pub(super) fn run(args: &[OsString], json_mode: bool) -> CommandOutput {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return CommandOutput {
            exit_code: exit::SUCCESS,
            stdout: if json_mode {
                format!(
                    "{{\"schema\":{},\"help\":{}}}\n",
                    quote(RESULT_SCHEMA),
                    quote(HELP)
                )
            } else {
                HELP.into()
            },
            stderr: String::new(),
        };
    }
    if args.len() != 1 {
        return diagnostic(exit::USAGE, &bad(HELP), json_mode);
    }
    let mut text = String::new();
    let read = File::open(std::path::Path::new(&args[0]))
        .and_then(|file| file.take(MAX_INPUT_BYTES + 1).read_to_string(&mut text));
    if let Err(error) = read {
        return diagnostic(exit::INPUT, &bad(error.to_string()), json_mode);
    }
    if text.len() as u64 > MAX_INPUT_BYTES {
        return diagnostic(exit::INPUT, &bad("scene exceeds 16 MiB"), json_mode);
    }
    let base = std::path::Path::new(&args[0])
        .parent()
        .map_or_else(std::path::PathBuf::new, std::path::Path::to_path_buf);
    let root = match J::parse(&text) {
        Ok(root) => root,
        Err(e) => {
            return diagnostic(
                exit::REFUSED,
                &bad(format!("invalid JSON: {e:?}")),
                json_mode,
            );
        }
    };
    let result = if root.get("study").is_some() {
        study::run(&root, &base, json_mode)
    } else {
        match Scene::from_root(&root, &base) {
            Ok(scene) => execute_within_budget(&scene, json_mode),
            Err(failure) => {
                let class = if failure.code == "cooling-cht-budget" {
                    exit::BUDGET
                } else {
                    exit::REFUSED
                };
                return diagnostic(class, &failure, json_mode);
            }
        }
    };
    match result {
        Ok(stdout) => CommandOutput {
            exit_code: exit::SUCCESS,
            stdout,
            stderr: String::new(),
        },
        Err(failure) => {
            let class = match failure.code {
                "cooling-cht-cancelled" => exit::BUDGET,
                "cooling-cht-input" => exit::REFUSED,
                "cooling-cht-not-steady" => exit::BUDGET,
                _ => exit::REFUSED,
            };
            diagnostic(class, &failure, json_mode)
        }
    }
}
