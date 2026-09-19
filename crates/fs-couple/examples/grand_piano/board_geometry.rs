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

use super::linear::BoardMode;
use fs_math::det;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateModel,
    PlateSection, SliceOptions, Stiffener,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

pub const HEADER: &str = "frankensim-board-geometry-si-v1";
// Matches board::read and the small dense bridge owner. Never truncate a
// certified slice silently just because the runtime budget is smaller.
const MAX_BOARD_MODES: usize = 32;
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
}

#[derive(Debug)]
pub struct PreparedBoard {
    pub modes: Vec<BoardMode>,
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
        })
    }

    /// Compute the complete certified slice (0, upper_hz]. Refuse an over-budget
    /// slice instead of silently dropping modes, or generating arbitrarily
    /// placed oscillators for missing geometry. All retained modes have SI,
    /// mass-normalized work-conjugate bridge force/displacement projections.
    pub fn prepare(&self, keys: &[u8], upper_hz: f64) -> Result<PreparedBoard, String> {
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
        let model = self.chart.assemble(&self.stiffeners, &self.supports).map_err(|e| e.to_string())?;
        if model.free == 0 { return Err("all soundboard degrees of freedom are constrained".into()); }
        let report = fs_plate::modes(&model, (0.0, (TAU * upper_hz).powi(2)), &SliceOptions::default())
            .map_err(|e| e.to_string())?;
        if report.below_low != 0 { return Err("unstable soundboard has negative stiffness eigenvalues".into()); }
        if report.modes.is_empty() || report.modes.len() > MAX_BOARD_MODES {
            return Err(format!("certified band contains {} modes; runtime admits 1..={MAX_BOARD_MODES}; choose a narrower explicit band", report.modes.len()));
        }
        let mesh = &self.chart.mesh;
        let mut volume_weights = vec![0.0; mesh.nodes.len()];
        for tri in &mesh.tris {
            let area = triangle_area(mesh, *tri);
            for &node in tri { volume_weights[node] += area / 3.0; }
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
                bridge[usize::from(site.midi - 21)] = (0..3).map(|i| site.weights[i] * w(tri[i])).sum();
            }
            let volume: f64 = volume_weights.iter().enumerate().map(|(node, area)| area * w(node)).sum();
            if !volume.is_finite() || bridge.iter().any(|x| !x.is_finite()) {
                return Err("soundboard port projection overflow".into());
            }
            modes.push(BoardMode {
                frequency_hz: det::sqrt(pair.lambda) / TAU,
                damping_ratio: self.damping_ratio, bridge, volume,
            });
            intervals.push((det::sqrt(pair.interval.0) / TAU, det::sqrt(pair.interval.1) / TAU));
        }
        let mass = self.mass_kg();
        if !mass.is_finite() || mass <= 0.0 { return Err("board mass overflow".into()); }
        Ok(PreparedBoard {
            modes, provenance: self.provenance.clone(), area_m2: mesh.total_area(),
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

fn nodal_displacement(model: &PlateModel, phi: &[f64], node: usize) -> f64 {
    model.dof_map[3 * node].map_or(0.0, |i| phi[i])
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
}
