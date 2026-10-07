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
    ChtError, EnergyConfig, EnergySolution, FanCurve, FanInlet, FluidProperties, FvBoundary,
    FvBuoyancyConfig, FvFlow, SimpleConfig, SolidMaterial, ThermalFace, ThermalSetup, Voxel,
    VoxelDomain, fv_natural_convection, simple_flow, solve_energy,
};
use fs_rep_mesh::{Soup, WindingOctree, winding_exact};
use json::JsonValue as J;

const MAX_INPUT_BYTES: u64 = 16 * 1024 * 1024;
const MAX_CELLS: usize = 4_000_000;
const MAX_STL_BYTES: u64 = 256 * 1024 * 1024;
const SCHEMA: &str = "frankensim.cooling-cht.v1";
const RESULT_SCHEMA: &str = "frankensim.cooling-cht.result.v1";
const NO_CLAIM: &str = "steady laminar constant-property flow on a staircase voxel grid at one declared resolution (no mesh-convergence claim); Boussinesq buoyancy only when gravity is declared; no turbulence model, radiation, or temperature-dependent properties; power-law convection is first order at high cell Peclet numbers; Estimated numerical evidence, not validated hardware or a ledger-backed .fsim run";
const HELP: &str = "Usage: frankensim [--json] cooling-cht <scene.json>\n\nSolve steady voxel conjugate heat transfer: finite-volume SIMPLEC airflow\n(forced, or natural/mixed with the Boussinesq force when gravity_m_s2 is\ndeclared) and one conservative energy equation over fluid and solid cells.\nThe scene declares size_m and voxel_m, a fluid (\"dry-air-300k\" or explicit\nproperties), materials, solids (boxes, or closed STL meshes placed by\nscale and offset_m; later solids override earlier ones),\nheat-source boxes (power spread over the solid cells they cover), and one\nrule per face x-, x+, y-, y+, z-, z+: inlet (velocity_m_s, temperature_k),\nfan (curve [[flow_m3_s, pressure_pa], ...], temperature_k; the flow is the\noperating point against the system), opening (ambient_k; pressure zero, flow either way), symmetry, or wall\n(adiabatic, or temperature_k, heat_flux_w_m2, or htc_w_m2_k with ambient_k).\nMissing faces are adiabatic walls. Request schema: frankensim.cooling-cht.v1.\nResults are Estimated single-resolution numerical evidence.\n";

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
    dims: [usize; 3],
    dx: f64,
    fluid: FluidProperties,
    materials: Vec<SolidMaterial>,
    solids: Vec<(u16, Shape)>,
    sources: Vec<Source>,
    faces: [FaceRule; 6],
    gravity: Option<[f64; 3]>,
    expansion: Option<f64>,
    reference: Option<f64>,
    tolerance: f64,
    max_iterations: usize,
    wall_seconds: f64,
}

const FACE_KEYS: [&str; 6] = ["x-", "x+", "y-", "y+", "z-", "z+"];

impl Scene {
    #[allow(clippy::too_many_lines)] // one schema, field by field
    fn parse(text: &str, base: &std::path::Path) -> Result<Self> {
        let root = J::parse(text).map_err(|e| bad(format!("invalid JSON: {e:?}")))?;
        if root.str_field("schema") != Some(SCHEMA) {
            return Err(bad(format!("schema must be {SCHEMA}")));
        }
        let dx = number(&root, "voxel_m", "scene")?;
        if dx <= 0.0 {
            return Err(bad("scene.voxel_m must be positive"));
        }
        let size = vec3(&root, "size_m", "scene")?;
        let mut dims = [0usize; 3];
        for a in 0..3 {
            let cells = size[a] / dx;
            let rounded = cells.round();
            if rounded < 1.0 || (cells - rounded).abs() > 1e-6 * rounded.max(1.0) {
                return Err(bad(format!(
                    "scene.size_m[{a}] = {} is not a whole number of {dx} m voxels",
                    size[a]
                )));
            }
            dims[a] = rounded as usize;
        }
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
            names.push(name);
            materials.push(SolidMaterial::new(
                name,
                number(item, "conductivity_w_m_k", &at)?,
            ));
        }
        if materials.len() > usize::from(u16::MAX) {
            return Err(bad("at most 65535 materials"));
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
            let shape = if item.get("stl").is_some() {
                Shape::Mesh(Box::new(MeshSolid::load(item, &at, base)?))
            } else {
                Shape::Box(Aabb::parse(item, &at)?)
            };
            solids.push((u16::try_from(index).expect("bounded above"), shape));
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
                Some("fan") => {
                    let points = field(rule, "curve", &at)?
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
                    FaceRule::Fan {
                        curve: FanCurve::new(&points)
                            .map_err(|e| bad(format!("{at}.curve: {e}")))?,
                        temperature: number(rule, "temperature_k", &at)?,
                    }
                }
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
        let wall_seconds = match root.get("limits") {
            Some(l) => optional_number(l, "wall_seconds", "limits")?.unwrap_or(3600.0),
            None => 3600.0,
        };
        if !(wall_seconds > 0.0) {
            return Err(bad("limits.wall_seconds must be positive"));
        }
        Ok(Self {
            dims,
            dx,
            fluid,
            materials,
            solids,
            sources,
            faces,
            gravity,
            expansion: optional_number(&root, "expansion_per_k", "scene")?,
            reference: optional_number(&root, "reference_temperature_k", "scene")?,
            tolerance,
            max_iterations: max_iterations as usize,
            wall_seconds,
        })
    }
}

#[allow(clippy::too_many_lines)] // build, solve, and one linear report
fn execute(scene: &Scene, gate: &CancelGate, json_mode: bool) -> Result<String> {
    let [nx, ny, nz] = scene.dims;
    let tie = 1e-6 * scene.dx;
    let domain = VoxelDomain::from_fn(nx, ny, nz, scene.dx, |p| {
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
                .map(|a| scene.dims[a] as f64 * scene.dx)
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
    let mut source_cells = Vec::with_capacity(scene.sources.len());
    for source in &scene.sources {
        let count = setup.add_uniform_power(&domain, source.power_w, |p| {
            source.region.contains(probe(p, tie))
                && !domain.is_fluid({
                    let cell = |v: f64| (v / scene.dx).floor() as usize;
                    domain.index(cell(p[0]), cell(p[1]), cell(p[2]))
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
    flow_config.max_iterations = scene.max_iterations;
    flow_config.fan = fans.first().map(|&side| FanInlet {
        face: Face3::ALL[side],
        curve: match scene.faces[side] {
            FaceRule::Fan { curve, .. } => curve,
            _ => unreachable!("filtered to fan faces"),
        },
    });
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
            let run = fv_natural_convection(
                &domain,
                &scene.fluid,
                &scene.materials,
                &setup,
                &config,
                gate,
            )
            .map_err(|e| solver_failure(&e))?;
            (run.flow, run.energy, Some(run.report.couplings))
        } else {
            let flow = simple_flow(&domain, &scene.fluid, &flow_config, gate)
                .map_err(|e| solver_failure(&e))?;
            let energy = solve_energy(
                &domain,
                &scene.fluid,
                &scene.materials,
                &flow.field,
                &setup,
                &EnergyConfig::default(),
                gate,
            )
            .map_err(|e| solver_failure(&e))?;
            (flow, energy, None)
        };
    let wall_s = started.elapsed().as_secs_f64();
    let temperature = &energy.temperature;
    // Per-material and per-source temperatures.
    let mut material_rows = Vec::new();
    for (index, material) in scene.materials.iter().enumerate() {
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
    let solver_name = if couplings.is_some() {
        "fv-simplec-boussinesq"
    } else {
        "fv-simplec"
    };
    if json_mode {
        let mut out = format!(
            "{{\"schema\":{},\"status\":\"completed\",\"solver\":{},\"cells\":{},\"fluid_cells\":{},\"voxel_m\":{}",
            quote(RESULT_SCHEMA),
            quote(solver_name),
            domain.cell_count(),
            domain.fluid_count(),
            num(scene.dx)?
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
        if let Some(couplings) = couplings {
            let _ = write!(out, ",\"energy_couplings\":{couplings}");
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
        let _ = write!(
            out,
            "}},\"energy\":{{\"iterations\":{},\"relative_residual\":{},\"source_w\":{},\"boundary_outflow_w\":{},\"advective_outflow_w\":{},\"balance_relative_residual\":{},\"max_cell_peclet\":{}}}",
            energy.report.iterations,
            num(energy.report.relative_residual)?,
            num(b.source_w)?,
            num(b.boundary_outflow_w)?,
            num(b.advective_outflow_w)?,
            num(b.relative_residual)?,
            num(energy.report.max_cell_peclet)?
        );
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
        if let Some((q, dp, _)) = r.fan {
            let _ = writeln!(out, "fan_flow_m3_s={q:e}\nfan_pressure_pa={dp}");
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
    let scene = match Scene::parse(&text, &base) {
        Ok(scene) => scene,
        Err(failure) => {
            let class = if failure.code == "cooling-cht-budget" {
                exit::BUDGET
            } else {
                exit::REFUSED
            };
            return diagnostic(class, &failure, json_mode);
        }
    };
    // The wall-time limit trips the same cancellation gate the solvers poll,
    // so an exhausted budget publishes nothing.
    let gate = CancelGate::new();
    let duration = Duration::from_secs_f64(scene.wall_seconds);
    let result = std::thread::scope(|scope| {
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
        let result = execute(&scene, &gate, json_mode);
        let _ = stop.send(());
        result
    });
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
