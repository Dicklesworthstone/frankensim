//! Paired upper/lower OBJ skins -> geometry-derived sections and crowned shell.
//! The explicit correspondence is a geometric hypothesis, not a thickness fit.
//! Existing OBJ, plate and shell owners still admit every output. No solver,
//! visual-material inference, mesh repair or external process is introduced.
use super::super::{Imported, Spec, MAX_OBJ_BYTES, MAX_SPEC_BYTES, MAX_NODES, MAX_TRIS,
    cross, dot, index, lower, number};
use super::super::super::crowned_board;
use fs_io::obj::{ObjDocument, ObjRegion, read_obj_document};
use std::collections::{BTreeMap, BTreeSet};

pub const HEADER: &str = "frankensim-board-skins-v1";

struct Skins {
    lower: String,
    pairs: BTreeMap<usize, usize>,
    thickness: [f64; 2],
}
impl Skins {
    fn read(text: &str) -> Result<Self, String> {
        if text.len() > MAX_SPEC_BYTES { return Err("skin correspondence exceeds 8 MiB".into()); }
        let mut rows = text.lines().map(str::trim).filter(|s| !s.is_empty() && !s.starts_with('#'));
        if rows.next() != Some(HEADER) { return Err(format!("expected {HEADER}")); }
        let mut s = Self { lower: String::new(), pairs: BTreeMap::new(), thickness: [0.; 2] };
        let mut seen = BTreeSet::new();
        let mut used_lower = BTreeSet::new();
        for row in rows {
            if row.len() > 8192 { return Err("skin row exceeds 8192 bytes".into()); }
            let f: Vec<_> = row.split(',').map(str::trim).collect();
            if f[0] != "pair" && !seen.insert(f[0].to_owned()) {
                return Err(format!("duplicate skin {} row", f[0]));
            }
            match f.as_slice() {
                ["lower", label] if !label.is_empty() => s.lower = (*label).into(),
                ["thickness", "geometry"] => {},
                ["thickness-range", lo, hi] => s.thickness = [number(lo)?, number(hi)?],
                ["pair", upper, lower] => {
                    let (a, b) = (index(upper)?, index(lower)?);
                    if s.pairs.len() >= MAX_NODES || s.pairs.insert(a, b).is_some()
                        || !used_lower.insert(b) {
                        return Err("skin pairs must be one-to-one and within the node budget".into());
                    }
                }
                _ => return Err(format!("unknown skin row or wrong field count: {}", f[0])),
            }
        }
        if s.lower.is_empty() || s.pairs.is_empty() || !seen.contains("thickness")
            || !(1e-5..=0.1).contains(&s.thickness[0])
            || !(s.thickness[0]..=0.1).contains(&s.thickness[1]) {
            return Err("require lower part, pairs, thickness,geometry and an SI thickness-range within 1e-5..0.1 m".into());
        }
        Ok(s)
    }
}

fn selected(doc: &ObjDocument, label: &str) -> Result<(Vec<usize>, BTreeSet<usize>), String> {
    let mut faces = Vec::new();
    let mut nodes = BTreeSet::new();
    for r in &doc.regions {
        if !r.has_label(label) { continue; }
        for face in r.triangles.clone() {
            if faces.len() >= MAX_TRIS { return Err("selected skin exceeds triangle budget".into()); }
            faces.push(face);
            for id in doc.soup.triangles[face] { nodes.insert(id as usize + 1); }
            if nodes.len() > MAX_NODES { return Err("selected skin exceeds node budget".into()); }
        }
    }
    if faces.is_empty() { return Err(format!("no faces in exact skin label {label:?}")); }
    Ok((faces, nodes))
}
fn chart(doc: &ObjDocument, id: usize, spec: &Spec) -> Result<[f64; 3], String> {
    let p = doc.soup.positions.get(id - 1).ok_or("skin vertex outside OBJ")?;
    let d = [(p.x-spec.origin[0])*spec.scale, (p.y-spec.origin[1])*spec.scale,
        (p.z-spec.origin[2])*spec.scale];
    let p = [dot(d, spec.u), dot(d, spec.v), dot(d, cross(spec.u, spec.v))];
    if p.iter().any(|x| !x.is_finite()) { return Err("skin coordinate transform overflow".into()); }
    Ok(p)
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { std::array::from_fn(|i| a[i]-b[i]) }
fn unit_triangle(p: [[f64; 3]; 3]) -> Result<[f64; 3], String> {
    let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let length = dot(n, n).sqrt();
    if !length.is_finite() || length <= 1e-16 { return Err("degenerate skin triangle".into()); }
    Ok(n.map(|x| x/length))
}

/// Upper part/material/support/beam identities come from the existing FSPI.
/// FSPS supplies the lower part and one-based upper->lower vertex pairs.
/// `thickness,geometry` explicitly replaces each material's nominal thickness
/// with mean normal separation at that triangle's three paired vertices.
/// Taper is piecewise constant per element, not an exact variable-section solid.
/// Side walls and unrelated cabinet geometry are deliberately not selected.
/// This does not prove global self-intersection freedom or specimen fidelity.
pub fn import(obj: &str, specification: &str, correspondence: &str) -> Result<Imported, String> {
    if obj.len() > MAX_OBJ_BYTES { return Err("OBJ exceeds 32 MiB".into()); }
    let mut spec = Spec::read(specification)?;
    let skins = Skins::read(correspondence)?;
    if skins.lower == spec.part { return Err("upper and lower skin labels must differ".into()); }
    let mut doc = read_obj_document(obj).map_err(|e| e.to_string())?;
    let (upper_faces, upper_nodes) = selected(&doc, &spec.part)?;
    let (lower_faces, lower_nodes) = selected(&doc, &skins.lower)?;
    if !upper_nodes.is_disjoint(&lower_nodes)
        || skins.pairs.keys().copied().collect::<BTreeSet<_>>() != upper_nodes
        || skins.pairs.values().copied().collect::<BTreeSet<_>>() != lower_nodes {
        return Err("pairs must cover exactly both disjoint selected skin vertex sets".into());
    }
    let inverse: BTreeMap<_, _> = skins.pairs.iter().map(|(&a, &b)| (b, a)).collect();
    let mut upper_topology = BTreeSet::new();
    let mut lower_topology = BTreeSet::new();
    for &face in &upper_faces {
        let mut t = doc.soup.triangles[face].map(|v| v as usize + 1);
        t.sort_unstable();
        if !upper_topology.insert(t) { return Err("duplicate upper skin face".into()); }
    }
    for &face in &lower_faces {
        let mut t = doc.soup.triangles[face].map(|v| inverse[&(v as usize + 1)]);
        t.sort_unstable();
        if !lower_topology.insert(t) { return Err("duplicate lower skin face".into()); }
    }
    if upper_topology != lower_topology {
        return Err("paired skins must have identical triangle connectivity; remesh or supply corresponding skins explicitly".into());
    }
    // Each entry is [midpoint, upper, lower], all in the declared SI chart.
    let mut points = BTreeMap::new();
    for (&a, &b) in &skins.pairs {
        let top = chart(&doc, a, &spec)?;
        let bottom = chart(&doc, b, &spec)?;
        let mid: [f64; 3] = std::array::from_fn(|i| 0.5*top[i]+0.5*bottom[i]);
        if mid[2].abs() > 0.05 { return Err("skin midpoint crown exceeds 50 mm".into()); }
        points.insert(a, [mid, top, bottom]);
    }
    let mut regions = Vec::with_capacity(upper_faces.len());
    let mut materials = BTreeMap::new();
    let mut minimum = f64::INFINITY;
    let mut maximum = 0.0_f64;
    let mut variation = 0.0_f64;
    for region in &doc.regions {
        if !region.has_label(&spec.part) { continue; }
        let name = region.material.as_deref().unwrap_or("");
        let card = spec.materials.get(name).or_else(|| spec.materials.get("*"))
            .ok_or_else(|| format!("upper skin has unmapped physical material {name:?}"))?;
        for face in region.triangles.clone() {
            let t = doc.soup.triangles[face].map(|v| v as usize + 1);
            let mid = t.map(|id| points[&id][0]);
            let raw = unit_triangle(mid)?;
            let sign = if raw[2] < 0.0 { -1.0 } else { 1.0 };
            let n = raw.map(|x| sign*x);
            if n[2] < 0.95 { return Err("paired midsurface is outside the shallow-shell slope limit".into()); }
            for side in [1, 2] {
                let direction = unit_triangle(t.map(|id| points[&id][side]))?;
                if sign*dot(direction, n) < 0.95 {
                    return Err("skin folds or diverges excessively from its midsurface".into());
                }
            }
            let mut h = [0.; 3];
            for i in 0..3 {
                let p = points[&t[i]];
                let gap = sub(p[1], p[2]);
                let length = dot(gap, gap).sqrt();
                h[i] = dot(gap, n);
                if !length.is_finite() || length <= 0.0 || h[i]/length < 0.95
                    || !h[i].is_finite() || h[i] < skins.thickness[0] || h[i] > skins.thickness[1] {
                    return Err(format!("skin pair at upper vertex {} is inverted, sheared or outside the declared normal-thickness range", t[i]));
                }
            }
            let mut section = *card;
            section[0] = h.iter().sum::<f64>()/3.0;
            let lo = h.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = h.iter().copied().fold(0., f64::max);
            minimum = minimum.min(lo); maximum = maximum.max(hi);
            variation = variation.max(hi-lo);
            let material = format!("solid_section_{face}");
            materials.insert(material.clone(), section);
            // Do not clone arbitrarily long source labels once per triangle.
            regions.push(ObjRegion { triangles: face..face+1,
                object: Some("solid_midsurface".into()), groups: Vec::new(), material: Some(material) });
        }
    }
    // Retain upper source identities for fixed nodes and bonded beam paths.
    // Flatten only the topology chart, then restore actual midpoint heights.
    for (&id, p) in &points {
        let mid = p[0];
        let plane: [f64; 3] = std::array::from_fn(|i|
            spec.origin[i]+(mid[0]*spec.u[i]+mid[1]*spec.v[i])/spec.scale);
        let v = &mut doc.soup.positions[id-1];
        v.x = plane[0]; v.y = plane[1]; v.z = plane[2];
    }
    doc.regions = regions;
    spec.materials = materials;
    spec.part = "solid_midsurface".into();
    let mut out = lower(&doc, &spec)?;
    let heights: Vec<_> = out.source_vertices.iter().map(|id| points[id][0][2]).collect();
    let description = format!("paired OBJ skins; actual midpoint crown and triangle-mean normal thickness; nominal FSPI thickness explicitly replaced; vertex normal thickness {minimum:.17e}..{maximum:.17e} m; maximum within-element thickness spread {variation:.17e} m; rectangular bonded beam sections unchanged");
    out.fsb = crowned_board::elevate(&out.fsb, &heights, &description)?;
    out.max_projection_m = heights.iter().map(|z| z.abs()).fold(0., f64::max);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    pub(super) const SPEC: &str = "frankensim-obj-board-v1\nsource,estimated,synthetic paired-skin regression\npart,upper\nunits,1\nframe,0,0,0,1,0,0,0,1,0\nflatness,1e-8\nmaterial,spruce,0.008,450,1e10,8e8,0.3,6e8,0.2\nsupport,clamped\nboundary,all\ndamping,0.01\npretension,0\nbridge,69,0.5,0.5\n";
    pub(super) const PAIRS: &str = "frankensim-board-skins-v1\nlower,lower\nthickness,geometry\nthickness-range,0.001,0.03\npair,1,6\npair,2,7\npair,3,8\npair,4,9\npair,5,10\n";
    pub(super) fn fixture(crown: f64, factor: f64) -> String {
        use std::fmt::Write;
        let mut obj = String::new();
        for side in [0.5, -0.5] {
            for (x,y,c) in [(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(0.5,0.5,crown)] {
                writeln!(obj,"v {x:.17e} {y:.17e} {:.17e}",c+side*factor*(0.006+0.003*x)).unwrap();
            }
        }
        obj.push_str("o upper\nusemtl spruce\nf 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\no lower\nusemtl ignored_shader\nf 6 10 7\nf 7 10 8\nf 8 10 9\nf 9 10 6\n");
        obj
    }
    #[test]
    fn skins_set_true_midpoint_crown_and_geometry_derived_structural_mass() {
        let obj = fixture(0.015,1.0);
        let a = import(&obj,SPEC,PAIRS).unwrap();
        let board = crowned_board::CrownedBoard::read(&a.fsb).unwrap();
        assert_eq!(a.source_vertices,vec![1,2,3,4,5]);
        assert!((board.max_height_m-0.015).abs()<1e-14);
        // Vertical paired columns: area * projected normal thickness is volume.
        assert!((board.mass_kg-450.0*0.0075).abs()<1e-10);
        assert_eq!(a.fsb,import(&obj,SPEC,PAIRS).unwrap().fsb);
        assert_eq!(a.fsb,import(&obj,&SPEC.replace("material,spruce,0.008,","material,spruce,0.020,"),PAIRS).unwrap().fsb);
    }
    #[test]
    fn imported_taper_changes_existing_shell_modes_not_just_metadata() {
        let a = import(&fixture(0.,1.),SPEC,PAIRS).unwrap();
        let b = import(&fixture(0.,2.),SPEC,PAIRS).unwrap();
        let a = crowned_board::CrownedBoard::read(&a.fsb).unwrap();
        let b = crowned_board::CrownedBoard::read(&b.fsb).unwrap();
        assert!((b.mass_kg-2.*a.mass_kg).abs()<1e-10);
        let ma = a.prepare(&[69],2000.).unwrap().modes;
        let mb = b.prepare(&[69],2000.).unwrap().modes;
        assert!(!ma.is_empty() && !mb.is_empty());
        assert!((ma[0].frequency_hz-mb[0].frequency_hz).abs()>1.0);
    }
    #[test]
    fn incomplete_inverted_or_ambiguous_pairs_never_fall_back_to_nominal_sections() {
        let obj = fixture(0.015,1.);
        for pairs in [PAIRS.replace("pair,5,10\n",""),format!("{PAIRS}pair,5,10\n"),
            PAIRS.replace("pair,5,10","pair,5,9"),PAIRS.replace("lower,lower","lower,upper"),
            PAIRS.replace("thickness,geometry\n",""),PAIRS.replace("0.001,0.03","0.012,0.03"),
            PAIRS.replace("pair,1,6","pair,1,7").replace("pair,2,7","pair,2,6")] {
            assert!(import(&obj,SPEC,&pairs).is_err());
        }
        assert!(import(&fixture(0.015,-1.),SPEC,PAIRS).is_err());
        assert!(import(&fixture(0.06,1.),SPEC,PAIRS).is_err());
        assert!(import(&obj.replace("f 6 10 7","f 6 8 7"),SPEC,PAIRS).is_err());
        assert!(import(&obj,&SPEC.replace("material,spruce","material,maple"),PAIRS).is_err());
        assert!(import(&obj,&SPEC.replace("pretension,0","pretension,10"),PAIRS).is_err());
    }
}
