//! Cold geometry -> fs-plate -> fs-modal -> reciprocal piano bridge ports.
//!
//! No modal frequency, modal mass, bridge shape, or radiating area is authored
//! here. They follow from the supplied mesh, sections, ribs and supports.
//! The admitted model is a FLAT, linearized, perfectly bonded thin plate:
//! crown, glue slip, rim compliance and anisotropic downbearing are not inferred.
//! Input must already be a conforming, non-overlapping mesh; this adapter
//! checks incidence and local quality, not general triangle intersection.
//! All expensive assembly, eigenanalysis and allocations precede rendering.
//!
//! This file format is intentionally not named "Steinway measurements". The
//! public Model D page publishes 9 -> 6 mm board thickness and wood species,
//! but not a complete outline, rib schedule, bridge map or elastic specimen
//! measurements: https://www.steinway.com/pianos/steinway/grand/model-d
//! The caller must supply those inputs and their authority explicitly.

use super::linear::{BoardMode, MAX_BOARD_MODES};
use fs_math::det;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateModel,
    PlateSection, SliceOptions, Stiffener,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

#[path = "board_motion.rs"]
pub mod motion;

pub const HEADER: &str = "frankensim-board-geometry-si-v1";
// Use the SAME retention budget as modal import, bridge mechanics and radiation.
// Never silently truncate the certified frequency slice to meet this budget.
const MAX_NODES: usize = 20_000;
const MAX_TRIANGLES: usize = 40_000;

#[derive(Clone, Debug)]
struct BridgeSite {
    midi: u8,
    triangle: usize,
    /// P1 contact patch at an explicitly supplied barycentric station.
    weights: [f64; 3],
}

#[derive(Debug)]
pub struct BoardGeometry {
    provenance: String,
    chart: PlateChart,
    stiffeners: Vec<Stiffener>,
    supports: AssemblyOptions,
    bridge_sites: Vec<BridgeSite>,
    /// Authored or measured modal damping ratio, not a fitted geometry result.
    damping_ratio: f64,
    /// Acoustic integration only; the structural mesh and eigensolve stay fixed.
    acoustic_refinement_levels: usize,
}

/// One degree-two triangle or subtriangle quadrature point for P1 displacement.
/// These are bare-board coordinates; the coupled bank owns their mass-loading
/// transformation. Areas sum to the panel area, not to its bounding rectangle.
#[derive(Clone, Debug)]
pub struct SurfaceSample {
    pub position_m: [f64; 3],
    pub area_m2: f64,
    pub mode_shape: Vec<f64>,
}

#[derive(Debug)]
pub struct PreparedBoard {
    pub surface: Vec<SurfaceSample>,
    pub modes: Vec<BoardMode>,
    /// Optional full-vector field from the same eigensolve. Ordinary Rayleigh
    /// preparation does not allocate it. No second modal basis is reconstructed.
    pub motion: Option<motion::MotionSurface>,
    pub provenance: String,
    pub area_m2: f64,
    /// Panel + stiffener physical mass, INCLUDING fixed boundary nodes.
    pub mass_kg: f64,
    /// Eigenvalue certificate converted to frequency intervals [Hz].
    pub frequency_intervals_hz: Vec<(f64, f64)>,
    pub free_dofs: usize,
}

fn scalar(fields: &[&str], i: usize) -> Result<f64, String> {
    let x = fields[i].parse::<f64>().map_err(|_| format!("invalid scalar field {}", i + 1))?;
    if !x.is_finite() { return Err(format!("nonfinite scalar field {}", i + 1)); }
    Ok(x)
}
fn index(fields: &[&str], i: usize) -> Result<usize, String> {
    fields[i].parse::<usize>().map_err(|_| format!("invalid index field {}", i + 1))
}
fn once<T>(slot: &mut Option<T>, value: T, tag: &str) -> Result<(), String> {
    if slot.is_some() { return Err(format!("duplicate {tag} row")); }
    *slot = Some(value);
    Ok(())
}

impl BoardGeometry {
    /// Strict tagged SI rows. Identifiers must be contiguous in input order.
    ///
    /// - source,measured|published|estimated|mixed,free text (commas allowed)
    /// - node,id,x_m,y_m
    /// - triangle,id,n0,n1,n2,h_m,rho_kg_m3,E_L_Pa,E_R_Pa,nu_LR,G_LR_Pa,grain_rad
    /// - support,clamped|simply_supported
    /// - fixed,node (explicit supports; no inferred rim constraint)
    /// - stiffener,E_Pa,G_Pa,A_m2,I_m4,J_m4,eccentricity_m,rho_kg_m3,node0,node1,...
    /// - bridge,midi,triangle_id,bary0,bary1,bary2
    /// - damping,dimensionless_ratio
    /// - pretension,N_per_m (uniform nonnegative membrane tension, NOT crown)
    ///
    /// Unknown rows, missing values, duplicate keys and nonphysical data refuse.
    pub fn read(text: &str) -> Result<Self, String> {
        let mut header = false;
        let mut source = None;
        let mut nodes = Vec::new();
        let mut triangles = Vec::new();
        let mut sections = Vec::new();
        let mut stiffeners = Vec::new();
        let mut sites = Vec::new();
        let mut fixed = BTreeSet::new();
        let mut keys = BTreeSet::new();
        let mut support = None;
        let mut damping = None;
        let mut pretension = None;
        for (line, raw) in text.lines().enumerate() {
            let raw = raw.trim();
            if raw.is_empty() || raw.starts_with('#') { continue; }
            if !header {
                if raw != HEADER { return Err(format!("line {}: expected {HEADER}", line + 1)); }
                header = true;
                continue;
            }
            let f: Vec<_> = raw.split(',').map(str::trim).collect();
            let row = (|| -> Result<(), String> {
                match (f[0], f.len()) {
                    ("source", n) if n >= 3 => {
                        if !["measured", "published", "estimated", "mixed"].contains(&f[1])
                            || f[2..].iter().all(|v| v.is_empty()) {
                            return Err("source needs an authority label and attribution".into());
                        }
                        once(&mut source, f[1..].join(","), "source")?;
                    }
                    ("node", 4) => {
                        if index(&f, 1)? != nodes.len() || nodes.len() >= MAX_NODES {
                            return Err("node ids must be contiguous and below the node budget".into());
                        }
                        nodes.push((scalar(&f, 2)?, scalar(&f, 3)?));
                    }
                    ("triangle", 12) => {
                        if index(&f, 1)? != triangles.len() || triangles.len() >= MAX_TRIANGLES {
                            return Err("triangle ids must be contiguous and below the element budget".into());
                        }
                        triangles.push([index(&f, 2)?, index(&f, 3)?, index(&f, 4)?]);
                        sections.push(PlateSection::orthotropic_plane_stress_at_angle(
                            scalar(&f, 7)?, scalar(&f, 8)?, scalar(&f, 9)?, scalar(&f, 10)?,
                            scalar(&f, 5)?, scalar(&f, 6)?, scalar(&f, 11)?,
                        ).map_err(|e| e.to_string())?);
                    }
                    ("support", 2) => {
                        let value = match f[1] {
                            "clamped" => EdgeSupport::Clamped,
                            "simply_supported" => EdgeSupport::SimplySupported,
                            _ => return Err("unknown support; specify clamped or simply_supported".into()),
                        };
                        once(&mut support, value, "support")?;
                    }
                    ("fixed", 2) => {
                        if !fixed.insert(index(&f, 1)?) { return Err("duplicate fixed node".into()); }
                    }
                    ("stiffener", n) if n >= 10 => {
                        let s = Stiffener {
                            e: scalar(&f, 1)?, g: scalar(&f, 2)?, area: scalar(&f, 3)?,
                            inertia: scalar(&f, 4)?, torsion: scalar(&f, 5)?,
                            eccentricity: scalar(&f, 6)?, density: scalar(&f, 7)?,
                            nodes: (8..n).map(|i| index(&f, i)).collect::<Result<_, _>>()?,
                        };
                        if [s.e, s.g, s.area, s.density].iter().any(|x| *x <= 0.0)
                            || s.inertia < 0.0 || s.torsion < 0.0 {
                            return Err("invalid stiffener constants".into());
                        }
                        stiffeners.push(s);
                    }
                    ("bridge", 6) => {
                        let key = index(&f, 1)?;
                        if !(21..=108).contains(&key) || !keys.insert(key) {
                            return Err("bridge key must be unique and within 21..108".into());
                        }
                        let weights = [scalar(&f, 3)?, scalar(&f, 4)?, scalar(&f, 5)?];
                        if weights.iter().any(|x| *x < 0.0 || *x > 1.0)
                            || (weights.iter().sum::<f64>() - 1.0).abs() > 1e-12 {
                            return Err("bridge barycentric weights must be in [0,1] and sum to one".into());
                        }
                        sites.push(BridgeSite { midi: key as u8, triangle: index(&f, 2)?, weights });
                    }
                    ("damping", 2) => {
                        let ratio = scalar(&f, 1)?;
                        if !(0.0..1.0).contains(&ratio) { return Err("damping ratio must be in [0,1)".into()); }
                        once(&mut damping, ratio, "damping")?;
                    }
                    ("pretension", 2) => {
                        let t = scalar(&f, 1)?;
                        if t < 0.0 { return Err("compressive buckling/crown needs a shell model, not this flat adapter".into()); }
                        once(&mut pretension, t, "pretension")?;
                    }
                    _ => return Err(format!("unknown row or wrong field count: {}", f[0])),
                }
                Ok(())
            })();
            row.map_err(|e| format!("line {}: {e}", line + 1))?;
        }
        if fixed.is_empty() { return Err("explicit supports are required".into()); }
        if sites.is_empty() { return Err("at least one bridge station is required".into()); }
        let mesh = PlateMesh::from_unstructured(nodes, triangles).map_err(|e| e.to_string())?;
        validate_topology(&mesh)?;
        let section = *sections.first().ok_or("missing panel sections")?;
        let chart = PlateChart::with_boundary_and_regions(
            mesh, section, fixed.into_iter().collect(), vec![],
        ).map_err(|e| e.to_string())?.with_element_sections(sections).map_err(|e| e.to_string())?;
        for site in &sites {
            if site.triangle >= chart.mesh.tris.len() { return Err("bridge references a missing triangle".into()); }
        }
        // Validate paths here rather than relying on the broader plate API's
        // less restrictive handling of repeated/negative-density beam inputs.
        let mut beams = BTreeSet::new();
        for s in &stiffeners {
            for pair in s.nodes.windows(2) {
                if pair[0] >= chart.mesh.nodes.len() || pair[1] >= chart.mesh.nodes.len() || pair[0] == pair[1] {
                    return Err("invalid stiffener node path".into());
                }
                let edge = (pair[0].min(pair[1]), pair[0].max(pair[1]));
                if !beams.insert(edge) { return Err("duplicate stiffener segment would double mass and stiffness".into()); }
            }
        }
        Ok(Self {
            provenance: source.ok_or("missing source authority/attribution")?, chart, stiffeners,
            supports: AssemblyOptions {
                support: support.ok_or("missing support type")?,
                pretension: pretension.ok_or("missing explicit pretension (use zero for none)")?,
            },
            bridge_sites: sites, damping_ratio: damping.ok_or("missing explicit damping ratio")?,
            acoustic_refinement_levels: 0,
        })
    }

    /// Uniformly subdivide each P1 radiating triangle into 4^levels cells,
    /// with the existing positive degree-two rule on every cell. This changes
    /// observation quadrature alone, not the structural mesh or modal basis.
    /// The unchanged shared receiver admits at most 120,000 surface samples.
    /// Cubic displacement and crowned boards require their own matching rule.
    pub fn with_acoustic_refinement(mut self, levels: usize) -> Result<Self, String> {
        if levels > 3 || self.chart.mesh.tris.len() * 3 * 4usize.pow(levels as u32) > 120_000 {
            return Err("acoustic refinement exceeds level 3 or the 120000-point receiver budget".into());
        }
        self.acoustic_refinement_levels = levels;
        Ok(self)
    }

    /// Compute the complete certified slice (0, upper_hz]. Refuse an over-budget
    /// slice instead of silently dropping modes, or generating arbitrarily
    /// placed oscillators for missing geometry. All retained modes have SI,
    /// mass-normalized work-conjugate bridge force/displacement projections.
    pub fn prepare(&self, keys: &[u8], upper_hz: f64) -> Result<PreparedBoard, String> {
        self.prepare_inner(keys, upper_hz, false, false, false, false)
    }
    /// Use the same physical pencil in mass-equilibrated coordinates. This
    /// is an opt-in numerical trial; it does not alter geometry or materials.
    pub fn prepare_mass_equilibrated(&self, keys: &[u8], upper_hz: f64) -> Result<PreparedBoard, String> {
        self.prepare_inner(keys, upper_hz, false, true, false, false)
    }
    /// Use exact P1 transverse panel inertia; slope and beam inertia remain
    /// lumped. This is an opt-in structural discretization trial.
    pub fn prepare_consistent_transverse_mass(
        &self,
        keys: &[u8],
        upper_hz: f64,
        mass_equilibrated: bool,
    ) -> Result<PreparedBoard, String> {
        self.prepare_inner(keys, upper_hz, false, mass_equilibrated, true, false)
    }
    /// Integrate the opt-in cubic edge-compatible transverse field and use
    /// its matching reciprocal bridge and radiation-sample shapes. The
    /// currently reconstructed board still requires independent convergence
    /// and measured acoustic validation before this can be a default.
    pub fn prepare_edge_cubic_transverse_mass(
        &self, keys: &[u8], upper_hz: f64, mass_equilibrated: bool,
    ) -> Result<PreparedBoard, String> {
        self.prepare_inner(keys, upper_hz, false, mass_equilibrated, false, true)
    }
    /// Retain full-vector nodal motion for a finite acoustic skin. The same
    /// eigensolve supplies mechanics and acoustics, including repeated modes.
    pub fn prepare_with_motion(&self, keys: &[u8], upper_hz: f64) -> Result<PreparedBoard, String> {
        self.prepare_inner(keys, upper_hz, true, false, false, false)
    }
    /// Retain the same full-vector skin motion while solving the unchanged
    /// geometric board pencil in mass-equilibrated coordinates.
    pub fn prepare_with_motion_mass_equilibrated(&self, keys: &[u8], upper_hz: f64) -> Result<PreparedBoard, String> {
        self.prepare_inner(keys, upper_hz, true, true, false, false)
    }
    fn prepare_inner(
        &self,
        keys: &[u8],
        upper_hz: f64,
        retain_motion: bool,
        mass_equilibrated: bool,
        consistent_transverse_mass: bool,
        edge_cubic_transverse_mass: bool,
    ) -> Result<PreparedBoard, String> {
        if consistent_transverse_mass && edge_cubic_transverse_mass {
            return Err("choose one transverse panel mass law".into());
        }
        if edge_cubic_transverse_mass && self.acoustic_refinement_levels != 0 {
            return Err("acoustic refinement currently requires P1 transverse displacement".into());
        }
        if !upper_hz.is_finite() || upper_hz <= 0.0 || upper_hz > 80_000.0 {
            return Err("invalid soundboard frequency ceiling".into());
        }
        let mut seen = BTreeSet::new();
        for &key in keys {
            if !(21..=108).contains(&key) || !seen.insert(key) {
                return Err("admitted keys must be unique piano MIDI keys".into());
            }
            if !self.bridge_sites.iter().any(|b| b.midi == key) {
                return Err(format!("geometry has no bridge station for key {key}"));
            }
        }
        if keys.is_empty() { return Err("empty key set".into()); }
        let model = if edge_cubic_transverse_mass {
            self.chart
                .assemble_edge_cubic_transverse_mass(&self.stiffeners, &self.supports)
        } else if consistent_transverse_mass {
            self.chart
                .assemble_consistent_transverse_mass(&self.stiffeners, &self.supports)
        } else {
            self.chart.assemble(&self.stiffeners, &self.supports)
        }
        .map_err(|e| e.to_string())?;
        if model.free == 0 { return Err("all soundboard degrees of freedom are constrained".into()); }
        let modal_options = SliceOptions {
            mass_diagonal_equilibration: mass_equilibrated,
            ..SliceOptions::default()
        };
        let report = fs_plate::modes(&model, (0.0, (TAU * upper_hz).powi(2)), &modal_options)
            .map_err(|e| e.to_string())?;
        if report.below_low != 0 { return Err("unstable soundboard has negative stiffness eigenvalues".into()); }
        if report.modes.is_empty() || report.modes.len() > MAX_BOARD_MODES {
            return Err(format!("certified band contains {} modes; runtime admits 1..={MAX_BOARD_MODES}; choose a narrower explicit band", report.modes.len()));
        }
        let mesh = &self.chart.mesh;
        let mut volume_weights = vec![0.0; mesh.nodes.len()];
        let mut cubic_volume_shape = vec![0.0; model.free];
        for tri in &mesh.tris {
            let area = triangle_area(mesh, *tri);
            if edge_cubic_transverse_mass {
                let x=tri.map(|node|mesh.nodes[node].0);
                let y=tri.map(|node|mesh.nodes[node].1);
                let mean=fs_plate::edge_cubic_transverse_mean_shape(&x,&y);
                for local in 0..9 {
                    if let Some(index)=model.dof_map[3*tri[local/3]+local%3] {
                        cubic_volume_shape[index]+=area*mean[local];
                    }
                }
            } else {
                for &node in tri { volume_weights[node] += area / 3.0; }
            }
        }
        let mut modes = Vec::with_capacity(report.modes.len());
        let mut intervals = Vec::with_capacity(report.modes.len());
        let mut m_phi = vec![0.0; model.free];
        for (mode_id, pair) in report.modes.iter().enumerate() {
            if pair.phi.len() != model.free || pair.phi.iter().any(|x| !x.is_finite())
                || !pair.lambda.is_finite() || pair.lambda <= 0.0
                || !pair.interval.0.is_finite() || pair.interval.0 <= 0.0
                || !pair.interval.1.is_finite() || pair.interval.1 < pair.interval.0 {
                return Err("invalid or nonpositive certified soundboard eigenpair".into());
            }
            model.m.spmv(&pair.phi, &mut m_phi);
            let norm = pair.phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum::<f64>();
            if !norm.is_finite() || (norm - 1.0).abs() > 1e-7 {
                return Err("soundboard eigenvector is not mass normalized".into());
            }
            for other in &report.modes[..mode_id] {
                let product = other.phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum::<f64>();
                if !product.is_finite() || product.abs() > 1e-7 {
                    return Err("soundboard modes are not mass orthogonal".into());
                }
            }
            let w = |node| nodal_displacement(&model, &pair.phi, node);
            let mut bridge = [0.0; 88];
            for site in &self.bridge_sites {
                let tri = mesh.tris[site.triangle];
                bridge[usize::from(site.midi - 21)] = if edge_cubic_transverse_mass {
                    let shape=cubic_triangle_shape(mesh,tri,site.weights);
                    local_shape_displacement(&model,&pair.phi,tri,&shape)
                } else {
                    (0..3).map(|i| site.weights[i] * w(tri[i])).sum()
                };
            }
            let volume: f64 = if edge_cubic_transverse_mass {
                cubic_volume_shape.iter().zip(&pair.phi).map(|(shape,value)|shape*value).sum()
            } else {
                volume_weights.iter().enumerate().map(|(node, area)| area * w(node)).sum()
            };
            if !volume.is_finite() || bridge.iter().any(|x| !x.is_finite()) {
                return Err("soundboard port projection overflow".into());
            }
            modes.push(BoardMode {
                frequency_hz: det::sqrt(pair.lambda) / TAU,
                damping_ratio: self.damping_ratio, bridge, volume,
            });
            intervals.push((det::sqrt(pair.interval.0) / TAU, det::sqrt(pair.interval.1) / TAU));
        }
        // The chosen symmetric cubic interior rule makes this positive
        // three-point rule reproduce its exact area mean; fs-plate proves
        // that identity for every local displacement DOF.
        let quadrature = acoustic_triangle_weights(self.acoustic_refinement_levels);
        let mut surface = Vec::with_capacity(quadrature.len() * mesh.tris.len());
        for tri in &mesh.tris {
            let area = triangle_area(mesh, *tri) / quadrature.len() as f64;
            for &weights in &quadrature {
                let mut position = [0.0; 3];
                for i in 0..3 {
                    position[0] += weights[i] * mesh.nodes[tri[i]].0;
                    position[1] += weights[i] * mesh.nodes[tri[i]].1;
                }
                let cubic_shape=edge_cubic_transverse_mass
                    .then(||cubic_triangle_shape(mesh,*tri,weights));
                let shape = report.modes.iter().map(|pair| {
                    if let Some(coefficients)=&cubic_shape {
                        local_shape_displacement(&model,&pair.phi,*tri,coefficients)
                    } else {
                        (0..3).map(|i| weights[i]
                            *nodal_displacement(&model,&pair.phi,tri[i])).sum()
                    }
                }).collect();
                surface.push(SurfaceSample { position_m: position, area_m2: area, mode_shape: shape });
            }
        }
        let motion = if retain_motion {
            let geometry=fs_plate::ShellMesh::new(mesh.nodes.iter().map(|&(x,y)|[x,y,0.]).collect(),mesh.tris.clone())
                .map_err(|e|e.to_string())?;
            let shapes=report.modes.iter().map(|pair| (0..mesh.nodes.len()).map(|node| {
                let at=|c:usize|model.dof_map[3*node+c].map_or(0.,|i|pair.phi[i]);
                // DKT coordinates are slopes, NOT physical axial rotations.
                [0.,0.,at(0),at(2),-at(1),0.]
            }).collect()).collect();
            Some(motion::MotionSurface::new(geometry,shapes)?)
        } else {None};
        let mass = self.mass_kg();
        if !mass.is_finite() || mass <= 0.0 { return Err("board mass overflow".into()); }
        let mut provenance = self.provenance.clone();
        if consistent_transverse_mass {
            provenance.push_str("; exact P1 transverse panel mass; lumped slope/beam mass");
        }
        if edge_cubic_transverse_mass {
            provenance.push_str("; opt-in cubic edge-compatible panel mass and bridge/surface fields; lumped rotary/beam inertia");
        }
        if mass_equilibrated {
            provenance.push_str("; mass-diagonal solver equilibration");
        }
        if self.acoustic_refinement_levels != 0 {
            provenance.push_str(&format!("; P1 acoustic-only refinement level {} ({} surface samples)",
                self.acoustic_refinement_levels, surface.len()));
        }
        Ok(PreparedBoard {
            modes, surface, motion, provenance, area_m2: mesh.total_area(),
            mass_kg: mass, frequency_intervals_hz: intervals, free_dofs: model.free,
        })
    }

    fn mass_kg(&self) -> f64 {
        let mesh = &self.chart.mesh;
        let mut mass: f64 = match self.chart.section_field() {
            fs_plate::PlateSectionField::Uniform(s) => mesh.total_area() * s.thickness * s.density,
            fs_plate::PlateSectionField::PerElement(sections) => mesh.tris.iter().zip(sections)
                .map(|(&tri, s)| triangle_area(mesh, tri) * s.thickness * s.density).sum(),
        };
        for s in &self.stiffeners {
            for pair in s.nodes.windows(2) {
                let (a, b) = (mesh.nodes[pair[0]], mesh.nodes[pair[1]]);
                mass += s.density * s.area * det::sqrt((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2));
            }
        }
        mass
    }
}

/// Equal-area barycentric subdivision. Level zero deliberately returns the
/// original coordinates directly, preserving the original floating-point path.
fn acoustic_triangle_weights(levels: usize) -> Vec<[f64; 3]> {
    let rule = [[2.0/3.0, 1.0/6.0, 1.0/6.0],
        [1.0/6.0, 2.0/3.0, 1.0/6.0], [1.0/6.0, 1.0/6.0, 2.0/3.0]];
    if levels == 0 { return rule.to_vec(); }
    let mut cells = vec![[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]];
    for _ in 0..levels {
        let mut children = Vec::with_capacity(4 * cells.len());
        for [a,b,c] in cells {
            let midpoint = |x: [f64; 3], y: [f64; 3]| std::array::from_fn(|i| 0.5*(x[i]+y[i]));
            let ab = midpoint(a,b); let bc = midpoint(b,c); let ca = midpoint(c,a);
            children.extend([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]);
        }
        cells = children;
    }
    cells.iter().flat_map(|cell| rule.map(|weights|
        std::array::from_fn(|i| (0..3).map(|j| weights[j]*cell[j][i]).sum())))
        .collect()
}

fn nodal_displacement(model: &PlateModel, phi: &[f64], node: usize) -> f64 {
    model.dof_map[3 * node].map_or(0.0, |i| phi[i])
}
fn cubic_triangle_shape(mesh: &PlateMesh, tri: [usize; 3], weights: [f64; 3]) -> [f64; 9] {
    let x=tri.map(|node|mesh.nodes[node].0);
    let y=tri.map(|node|mesh.nodes[node].1);
    fs_plate::edge_cubic_transverse_shape(&x,&y,weights)
}
fn local_shape_displacement(model: &PlateModel, phi: &[f64], tri: [usize; 3], shape: &[f64; 9]) -> f64 {
    (0..9).map(|local|model.dof_map[3*tri[local/3]+local%3]
        .map_or(0.0,|index|shape[local]*phi[index])).sum()
}
fn triangle_area(mesh: &PlateMesh, tri: [usize; 3]) -> f64 {
    let (a, b, c) = (mesh.nodes[tri[0]], mesh.nodes[tri[1]], mesh.nodes[tri[2]]);
    0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0))
}
fn validate_topology(mesh: &PlateMesh) -> Result<(), String> {
    let mut tris = BTreeSet::new();
    let mut edges: BTreeMap<(usize, usize), (usize, i32)> = BTreeMap::new();
    let mut used = vec![false; mesh.nodes.len()];
    for tri in &mesh.tris {
        let mut canonical = *tri;
        canonical.sort_unstable();
        if !tris.insert(canonical) { return Err("duplicate triangle would double panel mass/stiffness".into()); }
        for i in 0..3 {
            let (a, b) = (tri[i], tri[(i + 1) % 3]);
            used[a] = true;
            let edge = edges.entry((a.min(b), a.max(b))).or_insert((0, 0));
            edge.0 += 1;
            edge.1 += if a < b { 1 } else { -1 };
            if edge.0 > 2 { return Err("nonmanifold plate edge".into()); }
        }
    }
    if edges.values().any(|&(count, winding)| count == 2 && winding != 0) {
        return Err("neighboring triangles have inconsistent edge winding".into());
    }
    if used.contains(&false) { return Err("isolated mesh node would give singular mass".into()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> String {
        let mesh = PlateMesh::rectangle(0.7, 0.45, 3, 3);
        let mut text = format!("{HEADER}\nsource,estimated,analytic rectangular regression NOT Steinway\nsupport,simply_supported\npretension,0\ndamping,0.01\n");
        for (i, (x, y)) in mesh.nodes.iter().enumerate() {
            text.push_str(&format!("node,{i},{x},{y}\n"));
        }
        for (i, t) in mesh.tris.iter().enumerate() {
            text.push_str(&format!("triangle,{i},{},{},{},0.006,450,11000000000,800000000,0.3,600000000,0\n", t[0], t[1], t[2]));
        }
        for node in mesh.boundary_nodes() { text.push_str(&format!("fixed,{node}\n")); }
        text.push_str("bridge,69,8,0.2,0.3,0.5\n");
        text
    }

    #[test]
    fn measured_geometry_rows_do_not_silently_fill_missing_inputs() {
        let text = fixture();
        let g = BoardGeometry::read(&text).unwrap();
        assert_eq!(g.bridge_sites.len(), 1);
        assert!((g.mass_kg() - 0.7 * 0.45 * 0.006 * 450.0).abs() < 1e-12);
        for bad in [
            text.replace("pretension,0\n", ""),
            text.replace("0.006,450", "NaN,450"),
            text.replace("bridge,69,8,0.2,0.3,0.5", "bridge,69,8,0.2,0.3,0.6"),
            format!("{text}bridge,69,8,0.2,0.3,0.5\n"),
            format!("{text}unknown,42\n"),
            text.replace("source,estimated", "source,guaranteed-Steinway"),
        ] { assert!(BoardGeometry::read(&bad).is_err()); }
        assert!(g.prepare(&[60], 300.0).is_err());
    }

    #[test]
    fn ribs_add_their_actual_geometric_mass() {
        let text = fixture();
        let panel = BoardGeometry::read(&text).unwrap().mass_kg();
        let rib = "stiffener,10000000000,600000000,0.0003,0.00000001,0.00000002,0.01,500,4,5,6,7\n";
        let g = BoardGeometry::read(&format!("{text}{rib}")).unwrap();
        assert!((g.mass_kg() - panel - 500.0 * 0.0003 * 0.7).abs() < 1e-12);
        assert!(BoardGeometry::read(&format!("{text}{rib}{rib}")).is_err());
        assert!(BoardGeometry::read(&format!("{text}{}", rib.replace(",500,", ",-500,"))).is_err());
    }

    #[test]
    fn bridge_force_and_velocity_are_work_conjugate() {
        let g = BoardGeometry::read(&fixture()).unwrap();
        let model = g.chart.assemble(&g.stiffeners, &g.supports).unwrap();
        let site = &g.bridge_sites[0];
        let tri = g.chart.mesh.tris[site.triangle];
        let v: Vec<_> = (0..model.free).map(|i| (i as f64 + 0.2).sin()).collect();
        let point_v: f64 = (0..3).map(|i| site.weights[i] * nodal_displacement(&model, &v, tri[i])).sum();
        let force = 123.4;
        let mut generalized = vec![0.0; model.free];
        for (i, &node) in tri.iter().enumerate() {
            if let Some(dof) = model.dof_map[3 * node] { generalized[dof] += site.weights[i] * force; }
        }
        let modal_work: f64 = generalized.iter().zip(&v).map(|(f, v)| f * v).sum();
        assert!((modal_work - force * point_v).abs() < 1e-12);
    }

    #[test]
    fn actual_geometry_changes_eigenfrequencies_and_mass_normalized_ports() {
        let text = fixture();
        let a = BoardGeometry::read(&text).unwrap().prepare(&[69], 300.0).unwrap();
        // All material stiffnesses x4 => frequencies x2, mass/shape unchanged.
        let stiffer = text.replace("11000000000,800000000,0.3,600000000", "44000000000,3200000000,0.3,2400000000");
        let b = BoardGeometry::read(&stiffer).unwrap().prepare(&[69], 600.0).unwrap();
        assert!(!a.modes.is_empty());
        assert_eq!(a.modes.len(), b.modes.len());
        for (x, y) in a.modes.iter().zip(&b.modes) {
            assert!((y.frequency_hz / x.frequency_hz - 2.0).abs() < 1e-6);
            assert!((y.bridge[48].abs() - x.bridge[48].abs()).abs() < 1e-5);
        }
        assert!((a.mass_kg - b.mass_kg).abs() < 1e-12);
    }
    #[test]
    fn surface_quadrature_retains_area_and_signed_modal_volume() {
        let p=BoardGeometry::read(&fixture()).unwrap().prepare(&[69],300.0).unwrap();
        let area=p.surface.iter().map(|s|s.area_m2).sum::<f64>();
        assert!((area-p.area_m2).abs()<1e-12);
        for (i,m) in p.modes.iter().enumerate() {
            let volume=p.surface.iter().map(|s|s.area_m2*s.mode_shape[i]).sum::<f64>();
            assert!((volume-m.volume).abs()<1e-12);
        }
    }
    #[test]
    fn acoustic_refinement_preserves_mechanics_area_and_signed_volume() {
        let original=BoardGeometry::read(&fixture()).unwrap().prepare(&[69],300.0).unwrap();
        for level in 0..=3 {
            let refined=BoardGeometry::read(&fixture()).unwrap().with_acoustic_refinement(level)
                .unwrap().prepare(&[69],300.0).unwrap();
            assert_eq!(refined.surface.len(),original.surface.len()*4usize.pow(level as u32));
            assert_eq!(refined.frequency_intervals_hz,original.frequency_intervals_hz);
            assert_eq!(refined.mass_kg,original.mass_kg); assert_eq!(refined.free_dofs,original.free_dofs);
            assert!((refined.surface.iter().map(|p|p.area_m2).sum::<f64>()-original.area_m2).abs()<1e-12);
            for (i,(a,b)) in original.modes.iter().zip(&refined.modes).enumerate() {
                assert_eq!(a.frequency_hz,b.frequency_hz); assert_eq!(a.bridge,b.bridge);
                assert_eq!(a.volume,b.volume); assert_eq!(a.damping_ratio,b.damping_ratio);
                let integrated=refined.surface.iter().map(|p|p.area_m2*p.mode_shape[i]).sum::<f64>();
                assert!((integrated-b.volume).abs()<1e-12);
            }
            if level==0 {
                assert_eq!(refined.provenance,original.provenance);
                for (a,b) in original.surface.iter().zip(&refined.surface) {
                    assert_eq!(a.position_m,b.position_m); assert_eq!(a.area_m2,b.area_m2);
                    assert_eq!(a.mode_shape,b.mode_shape);
                }
            }
        }
        assert!(BoardGeometry::read(&fixture()).unwrap().with_acoustic_refinement(4).is_err());
        let refined=BoardGeometry::read(&fixture()).unwrap().with_acoustic_refinement(1).unwrap();
        assert!(refined.prepare_edge_cubic_transverse_mass(&[69],300.0,false).is_err());
    }
    #[test]
    fn refined_acoustic_rule_integrates_independent_barycentric_moments() {
        for level in 0..=3 {
            let points=acoustic_triangle_weights(level);
            let mean=|f:fn([f64;3])->f64| points.iter().map(|&p|f(p)).sum::<f64>()/points.len() as f64;
            assert!(points.iter().all(|p|p.iter().all(|&x|x>0.0)
                && (p.iter().sum::<f64>()-1.0).abs()<1e-14));
            assert!((mean(|p|p[0])-1.0/3.0).abs()<1e-14);
            assert!((mean(|p|p[0]*p[0])-1.0/6.0).abs()<1e-14);
            assert!((mean(|p|p[0]*p[1])-1.0/12.0).abs()<1e-14);
        }
    }
    #[test]
    fn acoustic_point_budget_refuses_before_structural_preparation() {
        let source=super::super::steinway_d::build(32).unwrap();
        assert!(BoardGeometry::read(&source.geometry).unwrap().with_acoustic_refinement(3).is_err());
    }
    #[test]
    fn retained_motion_uses_the_same_modal_basis_without_changing_ordinary_preparation() {
        let g=BoardGeometry::read(&fixture()).unwrap();
        let old=g.prepare(&[69],300.).unwrap();
        let full=g.prepare_with_motion(&[69],300.).unwrap();
        assert!(old.motion.is_none());
        let motion=full.motion.as_ref().unwrap();
        for (a,b) in old.modes.iter().zip(&full.modes) {
            assert_eq!(a.frequency_hz,b.frequency_hz);assert_eq!(a.bridge,b.bridge);
        }
        for point in &full.surface {
            let row=motion.normal_weights(point.position_m,[0.,0.,1.],0.).unwrap();
            for (a,b) in row.iter().zip(&point.mode_shape) {assert!((a-b).abs()<1e-11);}
        }
    }

    /// Research probe only: inspect the complete FE pencil above the runtime
    /// budget before considering a score/bridge-local reduction. A one-frequency
    /// ranking is deliberately not an admission rule for playable dynamics.
    #[test]
    #[ignore = "expensive Model D high-band eigenanalysis on the reviewed build host"]
    fn model_d_high_band_bridge_port_participation() {
        use fs_math::c64::C64;
        let divisions=std::env::var("FS_PIANO_PROBE_MESH_DIVISIONS")
            .map_or(Ok(24),|value|value.parse::<usize>()).unwrap();
        assert!([22,24,28,32,36,40,48,56,64].contains(&divisions));
        let row_gap=std::env::var("FS_PIANO_PROBE_ROW_GAP")
            .ok().map(|value|value.parse::<f64>().unwrap());
        assert!(row_gap.is_none() || divisions>32);
        let preset=if let Some(gap)=row_gap {super::super::steinway_d::build_probe_refined(divisions,gap)}
            else if divisions<=32 {super::super::steinway_d::build(divisions)}
            else {super::super::steinway_d::build_probe(divisions)}.unwrap();
        if let Ok(path)=std::env::var("FS_PIANO_PROBE_DUMP_GEOMETRY") {
            use std::io::Write;
            let mut file=std::fs::OpenOptions::new().write(true).create_new(true)
                .open(&path).expect("new research geometry path");
            file.write_all(preset.geometry.as_bytes()).expect("research geometry write");
            println!("MODEL_D_GEOMETRY path={path} bytes={}",preset.geometry.len());
            if std::env::var("FS_PIANO_PROBE_EXPORT_ONLY").as_deref()==Ok("1") {return;}
        }
        let geometry=BoardGeometry::read(&preset.geometry).unwrap();
        println!("MODEL_D_MESH mesh={divisions} row_gap={row_gap:?} nodes={} triangles={}",
            geometry.chart.mesh.nodes.len(),geometry.chart.mesh.tris.len());
        assert_eq!(geometry.bridge_sites.len(),88);
        assert!(geometry.bridge_sites.iter().all(|site|
            site.weights.iter().any(|weight|weight.abs()<1e-8)));
        let mut quality=geometry.chart.mesh.tris.iter().enumerate().map(|(index,triangle)| {
            let points=triangle.map(|node|geometry.chart.mesh.nodes[node]);
            let side=|a:(f64,f64),b:(f64,f64)|(a.0-b.0).powi(2)+(a.1-b.1).powi(2);
            let sum=side(points[0],points[1])+side(points[1],points[2])+side(points[2],points[0]);
            let twice_area=((points[1].0-points[0].0)*(points[2].1-points[0].1)
                -(points[1].1-points[0].1)*(points[2].0-points[0].0)).abs();
            (2.*det::sqrt(3.)*twice_area/sum,index)
        }).collect::<Vec<_>>();
        quality.sort_by(|a,b|a.0.total_cmp(&b.0));
        println!("MESH_QUALITY mesh={divisions} triangles={} min={:.9e} p01={:.9e} median={:.9e}",
            quality.len(),quality[0].0,quality[quality.len()/100].0,quality[quality.len()/2].0);
        for &(shape,index) in quality.iter().take(3) {
            let triangle=geometry.chart.mesh.tris[index];
            println!("WORST_TRIANGLE mesh={divisions} index={index} quality={shape:.9e} nodes={triangle:?} points={:?}",
                triangle.map(|node|geometry.chart.mesh.nodes[node]));
        }
        let consistent_mass=std::env::var("FS_PIANO_PROBE_CONSISTENT_MASS").as_deref()==Ok("1");
        let edge_cubic_mass=std::env::var("FS_PIANO_PROBE_EDGE_CUBIC_MASS").as_deref()==Ok("1");
        let edge_cubic_port=std::env::var("FS_PIANO_PROBE_EDGE_CUBIC_PORT").as_deref()==Ok("1");
        assert!(!(consistent_mass && edge_cubic_mass));
        assert!(!edge_cubic_port || edge_cubic_mass);
        let model=if edge_cubic_mass {
            geometry.chart.assemble_edge_cubic_transverse_mass(&geometry.stiffeners,&geometry.supports)
        } else if consistent_mass {
            geometry.chart.assemble_consistent_transverse_mass(&geometry.stiffeners,&geometry.supports)
        } else {
            geometry.chart.assemble(&geometry.stiffeners,&geometry.supports)
        }.unwrap();
        let report=fs_plate::modes(&model,(0.,(TAU*2200.).powi(2)),&SliceOptions {
            mass_diagonal_equilibration:true,..SliceOptions::default()
        }).unwrap();
        assert_eq!(report.below_low,0);
        assert_eq!(report.expected,report.modes.len());
        println!("MODEL_D_PENCIL mesh={divisions} mass={} modes={}",
            if edge_cubic_mass {"edge_cubic"} else if consistent_mass {"exact_p1"} else {"lumped"},
            report.modes.len());
        assert!(report.modes.len()>100);
        // The original cap probe is still binding through mesh 56. Mesh 64
        // has a measured 127-mode default slice and is a convergence probe,
        // not an expected over-cap refusal.
        if !consistent_mass && !edge_cubic_mass && divisions<=56 {
            assert!(report.modes.len()>MAX_BOARD_MODES);
        }
        let retained=report.modes.len().min(MAX_BOARD_MODES);
        let port=|phi:&[f64],key:u8| {
            let site=geometry.bridge_sites.iter().find(|site|site.midi==key).unwrap();
            let triangle=geometry.chart.mesh.tris[site.triangle];
            if edge_cubic_port {
                let active=site.weights.iter().enumerate()
                    .filter(|(_,weight)|**weight>1e-8).map(|(i,_)|i).collect::<Vec<_>>();
                if active.len()==2 {
                    let (i,j)=(active[0],active[1]);
                    let (ni,nj)=(triangle[i],triangle[j]);
                    let t=site.weights[j]/(site.weights[i]+site.weights[j]);
                    let (pi,pj)=(geometry.chart.mesh.nodes[ni],geometry.chart.mesh.nodes[nj]);
                    let slope=|node:usize| {
                        let dx=model.dof_map[3*node+1].map_or(0.,|index|phi[index]);
                        let dy=model.dof_map[3*node+2].map_or(0.,|index|phi[index]);
                        dx*(pj.0-pi.0)+dy*(pj.1-pi.1)
                    };
                    return (2.*t*t*t-3.*t*t+1.)*nodal_displacement(&model,phi,ni)
                        +(t*t*t-2.*t*t+t)*slope(ni)
                        +(-2.*t*t*t+3.*t*t)*nodal_displacement(&model,phi,nj)
                        +(t*t*t-t*t)*slope(nj);
                }
            }
            (0..3).map(|i|site.weights[i]*nodal_displacement(&model,phi,triangle[i])).sum::<f64>()
        };
        for (drive,receive) in [(60,60),(69,69),(84,84),(84,69),(69,60),(84,60)] {
            for sample in 0..=120 {
                let hz=1200.+sample as f64*1000./120.;
                let omega=TAU*hz;
                let mobility=report.modes.iter().fold(C64::ZERO,|total,pair| {
                    let natural=det::sqrt(pair.lambda);
                    total+C64::new(0.,-omega)
                        .scale(port(&pair.phi,drive)*port(&pair.phi,receive))
                        /C64::new(natural*natural-omega*omega,
                            -2.*geometry.damping_ratio*natural*omega)
                });
                println!("MOBILITY_GRID mesh={divisions} drive={drive} receive={receive} hz={hz:.9} re={:.17e} im={:.17e}",mobility.re,mobility.im);
            }
        }
        for (drive,receive,hz) in [(84,84,2110.),(84,69,2110.),(69,69,880.)] {
            let omega=TAU*hz;
            let mut terms:Vec<_>=report.modes.iter().enumerate().map(|(index,pair)| {
                let natural=det::sqrt(pair.lambda);
                let denominator=C64::new(natural*natural-omega*omega,
                    -2.*geometry.damping_ratio*natural*omega);
                let transfer=C64::new(0.,-omega)
                    .scale(port(&pair.phi,drive)*port(&pair.phi,receive))/denominator;
                (index,natural/TAU,transfer)
            }).collect();
            let full=terms.iter().fold(C64::ZERO,|total,term|total+term.2);
            let low=terms.iter().filter(|term|term.1<=1200.)
                .fold(C64::ZERO,|total,term|total+term.2);
            assert!(full.re.is_finite() && full.im.is_finite() && full.abs()>0.);
            if drive==receive {assert!(full.re>=0.);}
            terms.sort_by(|a,b|b.2.abs().total_cmp(&a.2.abs()).then(a.0.cmp(&b.0)));
            if drive==84 && receive==69 {
                for (index,natural,transfer) in terms.iter().take(12) {
                    println!("PORT_MODE mesh={divisions} index={index} hz={natural:.9} re={:.17e} im={:.17e} magnitude={:.17e}",
                        transfer.re,transfer.im,transfer.abs());
                }
                for (index,natural,transfer) in terms.iter().filter(|term|(90..=100).contains(&term.0)) {
                    println!("CLUSTER_MODE mesh={divisions} index={index} hz={natural:.9} re={:.17e} im={:.17e}",
                        transfer.re,transfer.im);
                }
            }
            let total_abs=terms.iter().map(|term|term.2.abs()).sum::<f64>();
            let count_99=terms.iter().scan(0.,|sum,term|{*sum+=term.2.abs();Some(*sum)})
                .position(|sum|sum>=0.99*total_abs).unwrap()+1;
            let omitted_bound=terms[retained..].iter().map(|term|term.2.abs()).sum::<f64>();
            let reduced=terms[..retained].iter().fold(C64::ZERO,|total,term|total+term.2);
            println!("mesh={divisions} full_modes={} drive={drive} receive={receive} hz={hz} Y_re={:.9e} Y_im={:.9e} |Y|={:.9e} low_1200_relative_error={:.6} top128_relative_error={:.6} omitted_absolute_bound_relative={:.6} terms_for_99pct_absolute={count_99}",
                report.modes.len(),full.re,full.im,full.abs(),(full-low).abs()/full.abs(),
                (full-reduced).abs()/full.abs(),omitted_bound/full.abs());
        }
        let mut source_bounds:Vec<_>=report.modes.iter().enumerate().map(|(index,pair)| {
            let drive=port(&pair.phi,84).abs();
            let maximum_receiver=(21..=108).map(|key|port(&pair.phi,key).abs())
                .fold(0.0_f64,f64::max);
            let natural=det::sqrt(pair.lambda);
            (index,drive*maximum_receiver/(2.*geometry.damping_ratio*natural))
        }).collect();
        source_bounds.sort_by(|a,b|b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        let total_bound=source_bounds.iter().map(|entry|entry.1).sum::<f64>();
        let omitted=source_bounds[retained..].iter().map(|entry|entry.1).sum::<f64>();
        let count_99=source_bounds.iter().scan(0.,|sum,entry|{*sum+=entry.1;Some(*sum)})
            .position(|sum|sum>=0.99*total_bound).unwrap()+1;
        println!("C6 drive to every bridge: top128 omitted global absolute mobility bound={omitted:.9e} m/(N s), bound fraction={:.6}, terms_for_99pct_bound={count_99}",
            omitted/total_bound);
        let selected: BTreeSet<_>=source_bounds[..retained].iter().map(|entry|entry.0).collect();
        for receive in [84,69] {
            let omega=TAU*2110.;
            let mut full=C64::ZERO;
            let mut reduced=C64::ZERO;
            for (index,pair) in report.modes.iter().enumerate() {
                let natural=det::sqrt(pair.lambda);
                let term=C64::new(0.,-omega)
                    .scale(port(&pair.phi,84)*port(&pair.phi,receive))
                    /C64::new(natural*natural-omega*omega,
                        -2.*geometry.damping_ratio*natural*omega);
                full=full+term;
                if selected.contains(&index) {reduced=reduced+term;}
            }
            println!("C6 source-bound top128 at 2110 Hz to key {receive}: relative complex error={:.6}, common all-receiver bound/full={:.6}",
                (full-reduced).abs()/full.abs(),omitted/full.abs());
        }
    }

    /// Compare the complete bridge projection of nearby high-band modes on
    /// two conforming meshes. Nearest frequency alone can misidentify a mode
    /// when a tight cluster changes order or rotates its shape.
    #[test]
    #[ignore = "two expensive exact-P1 Model D eigenanalyses on the reviewed build host"]
    fn model_d_high_band_bridge_mode_match() {
        let mut banks=Vec::new();
        let mut fields=Vec::new();
        for divisions in [40,48] {
            let preset=super::super::steinway_d::build_probe(divisions).unwrap();
            let geometry=BoardGeometry::read(&preset.geometry).unwrap();
            assert!(geometry.bridge_sites.iter().map(|site|site.midi).eq(21..=108));
            let model=geometry.chart.assemble_consistent_transverse_mass(
                &geometry.stiffeners,&geometry.supports).unwrap();
            let report=fs_plate::modes(&model,(0.,(TAU*2200.).powi(2)),&SliceOptions {
                mass_diagonal_equilibration:true,..SliceOptions::default()
            }).unwrap();
            assert_eq!(report.below_low,0);
            assert_eq!(report.expected,report.modes.len());
            assert!(report.modes.len()>100);
            let bank=report.modes.iter().map(|mode| {
                let bridge=geometry.bridge_sites.iter().map(|site| {
                    let triangle=geometry.chart.mesh.tris[site.triangle];
                    (0..3).map(|corner|site.weights[corner]
                        *nodal_displacement(&model,&mode.phi,triangle[corner]))
                        .sum::<f64>()
                }).collect::<Vec<_>>();
                let norm=det::sqrt(bridge.iter().map(|value|value*value).sum::<f64>());
                assert!(norm>0. && norm.is_finite());
                let c6=bridge[84-21];
                let a4=bridge[69-21];
                (det::sqrt(mode.lambda)/TAU,bridge,norm,c6,a4)
            }).collect::<Vec<_>>();
            banks.push(bank);
            fields.push((geometry,model,report.modes[95].phi.clone()));
        }
        let (coarse,fine)=(&banks[0],&banks[1]);
        for i in 90..=100 {
            let mut matches=(90..=100).map(|j| {
                let dot=coarse[i].1.iter().zip(&fine[j].1)
                    .map(|(left,right)|left*right).sum::<f64>();
                (j,(dot/(coarse[i].2*fine[j].2)).abs())
            }).collect::<Vec<_>>();
            matches.sort_by(|a,b|b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            for &(j,score) in matches.iter().take(3) {
                println!("BRIDGE_MODE_MATCH old={i} old_hz={:.9} old_c6={:.17e} old_a4={:.17e} new={j} new_hz={:.9} new_c6={:.17e} new_a4={:.17e} abs_cosine={score:.9}",
                    coarse[i].0,coarse[i].3,coarse[i].4,fine[j].0,fine[j].3,fine[j].4);
            }
        }
        for key in 79..=89 {
            let station=key-21;
            println!("C6_NEIGHBOR key={key} old_mode95={:.17e} new_mode95={:.17e}",
                coarse[95].1[station],fine[95].1[station]);
        }
        let position=|(geometry,_,_):&(BoardGeometry,PlateModel,Vec<f64>),key:u8| -> [f64;2] {
            let site=geometry.bridge_sites.iter().find(|site|site.midi==key).unwrap();
            let triangle=geometry.chart.mesh.tris[site.triangle];
            std::array::from_fn(|axis| (0..3).map(|corner| {
                let node=geometry.chart.mesh.nodes[triangle[corner]];
                site.weights[corner]*if axis==0 {node.0} else {node.1}
            }).sum::<f64>())
        };
        let (c6,next)=(position(&fields[0],84),position(&fields[0],85));
        for key in 79..=89 {
            let a=position(&fields[0],key);
            let b=position(&fields[1],key);
            assert!((a[0]-b[0]).abs()<1e-10 && (a[1]-b[1]).abs()<1e-10);
        }
        let delta=[next[0]-c6[0],next[1]-c6[1]];
        let length=det::sqrt(delta[0]*delta[0]+delta[1]*delta[1]);
        assert!(length>0.);
        let normal=[-delta[1]/length,delta[0]/length];
        let sign=if coarse[95].1.iter().zip(&fine[95].1)
            .map(|(a,b)|a*b).sum::<f64>()<0. {-1.} else {1.};
        for fraction in [0.,0.25,0.5,0.75,1.] {
            for offset in [-0.02,0.,0.02] {
                let point=[c6[0]+fraction*delta[0]+offset*normal[0],
                    c6[1]+fraction*delta[1]+offset*normal[1]];
                let values=fields.iter().enumerate().map(|(index,(geometry,model,phi))| {
                    let stencil=fs_plate::loading::PlatePointStencil::locate(
                        &geometry.chart.mesh,model,point,
                        fs_plate::loading::PlateLoadBudget {
                            max_nodes:20_000,max_triangles:40_000,
                        }).unwrap();
                    let zero=vec![0.;phi.len()];
                    let displacement=stencil.sample(phi,&zero).unwrap().displacement;
                    if index==0 {displacement/coarse[95].2}
                    else {sign*displacement/fine[95].2}
                }).collect::<Vec<_>>();
                println!("C6_LOCAL_FIELD fraction={fraction:.2} offset_m={offset:.3} x_m={:.9} y_m={:.9} old={:.17e} new={:.17e}",
                    point[0],point[1],values[0],values[1]);
            }
        }
    }

    /// Independent smooth-plate reference at the same modal rank as the C6
    /// instability. It separates DKT/mass discretization error from the
    /// irregular source outline, beams and bridge stations.
    #[test]
    #[ignore = "six expensive high-band analytic-plate eigenanalyses"]
    fn high_band_dkt_reference_plate() {
        let (a,b,e,nu,h,rho)=(1.7,1.0,11.0e9,0.30,0.008,380.0);
        let section=PlateSection::isotropic(e,nu,h,rho).unwrap();
        let d=section.d[0];
        let mut exact=(1..=30).flat_map(|m| (1..=30).map(move |n| {
            let k=std::f64::consts::PI.powi(2)
                *((m*m) as f64/(a*a)+(n*n) as f64/(b*b));
            k*det::sqrt(d/(rho*h))/TAU
        })).collect::<Vec<_>>();
        exact.sort_by(f64::total_cmp);
        for nx in [40,48] {
            let ny=(nx as f64/a).round() as usize;
            let chart=PlateChart::from_mesh(
                PlateMesh::rectangle(a,b,nx,ny),section.clone()).unwrap();
            let opts=AssemblyOptions {
                pretension:0.,support:EdgeSupport::SimplySupported,
            };
            for (mass,model) in [
                ("lumped",chart.assemble(&[],&opts).unwrap()),
                ("exact_p1",chart.assemble_consistent_transverse_mass(&[],&opts).unwrap()),
                ("edge_cubic",chart.assemble_edge_cubic_transverse_mass(&[],&opts).unwrap()),
            ] {
                let report=fs_plate::modes(&model,(0.,(TAU*2400.).powi(2)),&SliceOptions {
                    mass_diagonal_equilibration:true,..SliceOptions::default()
                }).unwrap();
                assert_eq!(report.below_low,0);
                assert_eq!(report.expected,report.modes.len());
                assert!(report.modes.len()>100);
                for rank in [20,50,95,100] {
                    let actual=det::sqrt(report.modes[rank].lambda)/TAU;
                    println!("ANALYTIC_DKT nx={nx} ny={ny} mass={mass} rank={rank} exact_hz={:.9} actual_hz={actual:.9} relative_error={:.9}",
                        exact[rank],(actual-exact[rank]).abs()/exact[rank]);
                }
            }
        }
    }

}
