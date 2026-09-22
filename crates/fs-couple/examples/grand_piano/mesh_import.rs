//! Cold OBJ-to-plate adapter. Render labels select geometry; explicitly supplied
//! sections supply physics. No material is inferred from a shader or a filename.
//! The output goes through the SAME BoardGeometry admission used by grand_piano.
use super::board_geometry::BoardGeometry;
use fs_io::obj::{ObjDocument, read_obj_document};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

pub const HEADER: &str = "frankensim-obj-board-v1";
pub const MAX_OBJ_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_SPEC_BYTES: usize = 8 * 1024 * 1024;
const MAX_NODES: usize = 20_000;
const MAX_TRIS: usize = 40_000;

fn number(value: &str) -> Result<f64, String> {
    let n: f64 = value.parse().map_err(|_| format!("invalid scalar {value:?}"))?;
    if !n.is_finite() { return Err("nonfinite scalar".into()); }
    Ok(n)
}
fn index(value: &str) -> Result<usize, String> {
    value.parse::<usize>().ok().filter(|n| *n > 0)
        .ok_or_else(|| "OBJ vertex references must be positive, one-based indices".into())
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { (0..3).map(|i| a[i] * b[i]).sum() }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn area2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}

#[derive(Debug)]
struct Spec {
    source: String,
    part: String,
    scale: f64,
    origin: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
    flatness: f64,
    // h,rho,E_L,E_R,nu_LR,G_LR,grain_rad, in board-chart coordinates.
    materials: BTreeMap<String, [f64; 7]>,
    fixed: BTreeSet<usize>,
    boundary: bool,
    bridges: BTreeMap<u8, [f64; 2]>,
    // Existing FSB scalar columns and SOURCE OBJ vertex indices.
    stiffeners: Vec<([f64; 7], Vec<usize>)>,
    controls: Vec<String>,
}
impl Spec {
    fn read(text: &str) -> Result<Self, String> {
        if text.len() > MAX_SPEC_BYTES { return Err("board import specification exceeds 8 MiB".into()); }
        let mut rows = text.lines().enumerate().filter_map(|(n, row)| {
            let row = row.trim();
            (!row.is_empty() && !row.starts_with('#')).then_some((n + 1, row))
        });
        if rows.next().map(|(_, row)| row) != Some(HEADER) {
            return Err(format!("expected {HEADER}"));
        }
        let mut s = Self { source: String::new(), part: String::new(), scale: 0.0,
            origin: [0.0;3], u:[0.0;3], v:[0.0;3], flatness:0.0,
            materials:BTreeMap::new(), fixed:BTreeSet::new(), boundary:false,
            bridges:BTreeMap::new(), stiffeners:Vec::new(), controls:Vec::new() };
        let mut seen = BTreeSet::new();
        for (line, row) in rows {
            // Bound a row before collecting fields (beam paths may be long).
            if row.len() > 512 * 1024 { return Err(format!("line {line}: row exceeds budget")); }
            let f: Vec<_> = row.split(',').map(str::trim).collect();
            if f.len() > MAX_NODES + 8 { return Err(format!("line {line}: too many fields")); }
            let parsed = (|| -> Result<(), String> {
                if ["source","part","units","frame","flatness","support","damping","pretension","boundary"]
                    .contains(&f[0]) && !seen.insert(f[0].to_owned()) {
                    return Err(format!("duplicate {} row", f[0]));
                }
                match (f[0], f.len()) {
                    ("source", n) if n >= 3 => { s.source = f[1..].join(","); }
                    ("part", 2) if !f[1].is_empty() => { s.part = f[1].into(); }
                    ("units", 2) => { s.scale = number(f[1])?; }
                    ("flatness", 2) => { s.flatness = number(f[1])?; }
                    ("frame", 10) => {
                        for i in 0..3 {
                            s.origin[i] = number(f[i+1])?;
                            s.u[i] = number(f[i+4])?;
                            s.v[i] = number(f[i+7])?;
                        }
                    }
                    ("material", 9) if !f[1].is_empty() => {
                        let mut values = [0.0;7];
                        for (i, value) in values.iter_mut().enumerate() { *value = number(f[i+2])?; }
                        // The existing orthotropic owner, not a second D-matrix.
                        fs_plate::PlateSection::orthotropic_plane_stress_at_angle(
                            values[2], values[3], values[4], values[5], values[0], values[1], values[6])
                            .map_err(|e| e.to_string())?;
                        if s.materials.len() >= MAX_TRIS || s.materials.insert(f[1].into(), values).is_some() {
                            return Err("duplicate material or material budget exceeded".into());
                        }
                    }
                    ("fixed", 2) => {
                        if s.fixed.len() >= MAX_NODES || !s.fixed.insert(index(f[1])?) {
                            return Err("duplicate fixed vertex or support budget exceeded".into());
                        }
                    }
                    ("boundary", 2) if f[1] == "all" => { s.boundary = true; }
                    ("bridge", 4) => {
                        let key: u8 = f[1].parse().map_err(|_| "invalid bridge key")?;
                        if !(21..=108).contains(&key)
                            || s.bridges.insert(key, [number(f[2])?,number(f[3])?]).is_some() {
                            return Err("bridge keys must be unique and in 21..108".into());
                        }
                    }
                    ("stiffener", n) if n >= 10 => {
                        if s.stiffeners.len() >= MAX_TRIS { return Err("stiffener budget exceeded".into()); }
                        let mut values = [0.0;7];
                        for (i, value) in values.iter_mut().enumerate() { *value = number(f[i+1])?; }
                        let vertices = f[8..].iter().map(|v| index(v)).collect::<Result<Vec<_>,_>>()?;
                        s.stiffeners.push((values, vertices));
                    }
                    ("support" | "damping" | "pretension", 2) => s.controls.push(f.join(",")),
                    _ => return Err(format!("unknown row or wrong field count: {}", f[0])),
                }
                Ok(())
            })();
            parsed.map_err(|e| format!("line {line}: {e}"))?;
        }
        for required in ["source","part","units","frame","flatness","support","damping","pretension"] {
            if !seen.contains(required) { return Err(format!("missing {required} row")); }
        }
        if !(1e-9..=1e6).contains(&s.scale) || !(1e-12..=1e-3).contains(&s.flatness)
            || (dot(s.u,s.u)-1.0).abs() > 1e-10 || (dot(s.v,s.v)-1.0).abs() > 1e-10
            || dot(s.u,s.v).abs() > 1e-10 {
            return Err("units must be positive; frame axes orthonormal; flatness must be 1e-12..1e-3 metres".into());
        }
        if s.materials.is_empty() || s.bridges.is_empty() || (!s.boundary && s.fixed.is_empty()) {
            return Err("explicit physical materials, bridge stations and supports are required".into());
        }
        Ok(s)
    }
}

pub struct Imported {
    pub fsb: String,
    pub source_vertices: Vec<usize>,
    pub triangles: usize,
    pub max_projection_m: f64,
}

/// Leaves the imported surface Estimate-only. Local plate/topology admission
/// is not a global intersection certificate or a measured-instrument claim.
pub fn import(obj: &str, specification: &str) -> Result<Imported, String> {
    if obj.len() > MAX_OBJ_BYTES { return Err("OBJ exceeds 32 MiB".into()); }
    let spec = Spec::read(specification)?;
    let document = read_obj_document(obj).map_err(|e| e.to_string())?;
    lower(&document, &spec)
}
fn lower(doc: &ObjDocument, spec: &Spec) -> Result<Imported, String> {
    let mut selected = Vec::new();
    let mut source_vertices = BTreeSet::new();
    for region in &doc.regions {
        if !region.has_label(&spec.part) { continue; }
        let material = region.material.as_deref().unwrap_or("");
        let section = spec.materials.get(material).or_else(|| spec.materials.get("*"))
            .ok_or_else(|| format!("part {:?} has unmapped physical material {material:?}", spec.part))?;
        for face in region.triangles.clone() {
            if selected.len() == MAX_TRIS { return Err("selected panel exceeds 40000 triangles".into()); }
            let tri = doc.soup.triangles[face];
            for v in tri { source_vertices.insert(v as usize + 1); }
            if source_vertices.len() > MAX_NODES { return Err("selected panel exceeds 20000 vertices".into()); }
            selected.push((tri, *section));
        }
    }
    if selected.is_empty() { return Err(format!("no triangles in exact OBJ part {:?}", spec.part)); }
    let source_vertices: Vec<_> = source_vertices.into_iter().collect();
    let remap: BTreeMap<_,_> = source_vertices.iter().enumerate().map(|(i,&v)| (v,i)).collect();
    let normal = cross(spec.u,spec.v);
    let mut nodes = Vec::with_capacity(source_vertices.len());
    let mut occupied = BTreeSet::new();
    let mut max_projection_m = 0.0_f64;
    for &id in &source_vertices {
        let p = doc.soup.positions[id-1];
        let d = [(p.x-spec.origin[0])*spec.scale, (p.y-spec.origin[1])*spec.scale, (p.z-spec.origin[2])*spec.scale];
        let xy = [dot(d,spec.u),dot(d,spec.v)];
        let distance = dot(d,normal).abs();
        if xy.iter().any(|v| !v.is_finite()) || !distance.is_finite() || distance > spec.flatness {
            return Err(format!("OBJ vertex {id} lies outside the finite flat-panel chart; crown/solid geometry requires a shell"));
        }
        let bits = xy.map(|v| if v == 0.0 {0} else {v.to_bits()});
        if !occupied.insert(bits) { return Err("coincident projected vertices: weld/repair explicitly before import".into()); }
        max_projection_m = max_projection_m.max(distance);
        nodes.push(xy);
    }
    let mut out = format!("frankensim-board-geometry-si-v1\nsource,{}\n# OBJ selected part: {:?}; maximum plane projection {:.17e} m\n",
        spec.source, spec.part, max_projection_m);
    for (i,p) in nodes.iter().enumerate() { writeln!(out,"node,{i},{:.17e},{:.17e}",p[0],p[1]).unwrap(); }
    let mut triangles = Vec::new();
    let mut unique = BTreeSet::new();
    let mut edges = BTreeMap::<(usize,usize),usize>::new();
    for (i,(source,section)) in selected.iter().enumerate() {
        let mut t = source.map(|v| remap[&(v as usize + 1)]);
        let area = area2(nodes[t[0]],nodes[t[1]],nodes[t[2]]);
        if !area.is_finite() || area.abs() <= 1e-16 { return Err("degenerate or overflowing projected triangle".into()); }
        if area < 0.0 { t.swap(1,2); }
        let mut canonical = t; canonical.sort_unstable();
        if !unique.insert(canonical) { return Err("duplicate projected face".into()); }
        for (a,b) in [(t[0],t[1]),(t[1],t[2]),(t[2],t[0])] {
            *edges.entry((a.min(b),a.max(b))).or_default() += 1;
        }
        write!(out,"triangle,{i},{},{},{}",t[0],t[1],t[2]).unwrap();
        for value in section { write!(out,",{value:.17e}").unwrap(); } out.push('\n');
        triangles.push(t);
    }
    if edges.values().any(|&n| n > 2) { return Err("nonmanifold panel edge".into()); }
    let mut fixed = spec.fixed.iter().map(|v| remap.get(v).copied()
        .ok_or_else(|| format!("fixed OBJ vertex {v} is absent from selected panel")))
        .collect::<Result<BTreeSet<_>,_>>()?;
    if spec.boundary {
        for (&(a,b),&count) in &edges { if count == 1 { fixed.insert(a); fixed.insert(b); } }
    }
    for node in fixed { writeln!(out,"fixed,{node}").unwrap(); }
    for control in &spec.controls { writeln!(out,"{control}").unwrap(); }
    for (values,vertices) in &spec.stiffeners {
        out.push_str("stiffener");
        for value in values { write!(out,",{value:.17e}").unwrap(); }
        for id in vertices {
            let node = remap.get(id).ok_or_else(|| format!("stiffener OBJ vertex {id} is absent from panel"))?;
            write!(out,",{node}").unwrap();
        }
        out.push('\n');
    }
    for (&key,&p) in &spec.bridges {
        let (triangle,weights) = locate(p,&nodes,&triangles)
            .ok_or_else(|| format!("bridge key {key} is outside selected panel; nearest-node snapping is not allowed"))?;
        writeln!(out,"bridge,{key},{triangle},{:.17e},{:.17e},{:.17e}",weights[0],weights[1],weights[2]).unwrap();
    }
    // Real consumer admission, including orthotropy, supports, incidence,
    // bridge partition of unity and duplicate stiffener segments.
    BoardGeometry::read(&out)?;
    Ok(Imported { fsb:out, source_vertices, triangles:triangles.len(), max_projection_m })
}
fn locate(p:[f64;2],nodes:&[[f64;2]],triangles:&[[usize;3]]) -> Option<(usize,[f64;3])> {
    for (i,&[a,b,c]) in triangles.iter().enumerate() {
        let twice = area2(nodes[a],nodes[b],nodes[c]);
        let mut w = [area2(p,nodes[b],nodes[c])/twice,area2(nodes[a],p,nodes[c])/twice,area2(nodes[a],nodes[b],p)/twice];
        if w.iter().all(|v| v.is_finite() && (-1e-10..=1.0+1e-10).contains(v)) {
            for v in &mut w { *v = v.clamp(0.0,1.0); }
            let sum:f64 = w.iter().sum(); for v in &mut w { *v /= sum; }
            return Some((i,w));
        }
    }
    None
}

/// Export an admitted native FSB as an editable material-labelled midsurface
/// and import specification. Preserves every element section, explicit support,
/// beam section/path and bridge location; it does NOT emit a cabinet render mesh.
pub fn export(text:&str) -> Result<(String,String),String> {
    if text.len() > MAX_SPEC_BYTES { return Err("native board exceeds 8 MiB".into()); }
    BoardGeometry::read(text)?;
    let mut obj = String::from("# Physical midsurface in metres, not a measured digital twin\no soundboard\n");
    let mut spec = format!("{HEADER}\npart,soundboard\nunits,1\nframe,0,0,0,1,0,0,0,1,0\nflatness,1e-8\n");
    let rows:Vec<Vec<&str>> = text.lines().map(str::trim).filter(|r| !r.is_empty() && !r.starts_with('#'))
        .map(|r| r.split(',').map(str::trim).collect()).collect();
    let mut nodes = Vec::new(); let mut triangles = Vec::new();
    let mut materials = BTreeMap::<String,usize>::new();
    // OBJ positive indices reference vertices already declared; native FSB
    // permits interleaved record types, so emit all vertices first.
    for f in rows.iter().filter(|f| f[0] == "node") {
        let p = [number(f[2])?,number(f[3])?]; nodes.push(p);
        writeln!(obj,"v {:.17e} {:.17e} 0",p[0],p[1]).unwrap();
    }
    for f in rows.iter().filter(|f| f[0] == "triangle") {
                let t = [f[2].parse::<usize>().map_err(|_|"invalid triangle")?,
                    f[3].parse::<usize>().map_err(|_|"invalid triangle")?, f[4].parse::<usize>().map_err(|_|"invalid triangle")?];
                let signature = f[5..].join(","); let next = materials.len();
                let material = *materials.entry(signature.clone()).or_insert(next);
                // A material name references all seven physical section fields.
                writeln!(obj,"usemtl section_{material}\nf {} {} {}",t[0]+1,t[1]+1,t[2]+1).unwrap();
                triangles.push(t);
    }
    for (values,id) in materials { writeln!(spec,"material,section_{id},{values}").unwrap(); }
    for f in &rows {
        match f[0] {
            "source" | "support" | "damping" | "pretension" => { writeln!(spec,"{}",f.join(",")).unwrap(); }
            "fixed" => { writeln!(spec,"fixed,{}",f[1].parse::<usize>().map_err(|_|"invalid support")?+1).unwrap(); }
            "stiffener" => {
                spec.push_str(&f[..8].join(","));
                for id in &f[8..] { write!(spec,",{}",id.parse::<usize>().map_err(|_|"invalid beam node")?+1).unwrap(); }
                spec.push('\n');
            }
            "bridge" => {
                let id:usize = f[2].parse().map_err(|_|"invalid bridge triangle")?;
                let t = triangles[id]; let w = [number(f[3])?,number(f[4])?,number(f[5])?];
                let x:f64 = (0..3).map(|j|w[j]*nodes[t[j]][0]).sum();
                let y:f64 = (0..3).map(|j|w[j]*nodes[t[j]][1]).sum();
                writeln!(spec,"bridge,{},{x:.17e},{y:.17e}",f[1]).unwrap();
            }
            _ => {}
        }
    }
    // Enforce our interchange subset before publishing a pair that cannot
    // come back (e.g. duplicate planar vertices, overlong source labels).
    import(&obj,&spec)?;
    Ok((obj,spec))
}

#[cfg(test)]
mod tests {
    use super::*;
    const OBJ:&str="o decoration\nv 90 90 90\no soundboard\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0.5 0.5 0\nusemtl spruce\nf 2 3 6\nf 3 4 6\nf 4 5 6\nf 5 2 6\n";
    const SPEC:&str="frankensim-obj-board-v1\nsource,estimated,synthetic material/mesh regression\npart,soundboard\nunits,1\nframe,0,0,0,1,0,0,0,1,0\nflatness,1e-8\nmaterial,spruce,0.008,450,1e10,8e8,0.3,6e8,0.2\nsupport,simply_supported\nboundary,all\ndamping,0.01\npretension,0\nbridge,69,0.5,0.5\n";
    #[test]
    fn selected_material_mesh_reaches_the_real_modal_solver() {
        let imported = import(OBJ,SPEC).unwrap();
        assert_eq!(imported.source_vertices,vec![2,3,4,5,6]);
        assert_eq!(imported.triangles,4);
        let prepared = BoardGeometry::read(&imported.fsb).unwrap().prepare(&[69],1000.0).unwrap();
        assert!((prepared.area_m2-1.0).abs()<1e-12);
        assert!((prepared.mass_kg-3.6).abs()<1e-10);
        assert!(!prepared.modes.is_empty());
        assert_eq!(imported.fsb,import(OBJ,SPEC).unwrap().fsb);
    }
    #[test]
    fn explicit_frame_and_millimetres_are_covariant() {
        let rotated="o soundboard\nv 10 20 30\nv 10 1020 30\nv 10 1020 1030\nv 10 20 1030\nv 10 520 530\nusemtl spruce\nf 1 2 5\nf 2 3 5\nf 3 4 5\nf 4 1 5\n";
        let spec=SPEC.replace("units,1\n","units,0.001\n")
            .replace("frame,0,0,0,1,0,0,0,1,0","frame,10,20,30,0,1,0,0,0,1");
        assert_eq!(import(OBJ,SPEC).unwrap().fsb,import(rotated,&spec).unwrap().fsb);
    }
    #[test]
    fn malformed_inputs_do_not_turn_into_an_estimated_fallback() {
        for spec in [SPEC.replace("part,soundboard","part,missing"),
            SPEC.replace("material,spruce","material,maple"),
            SPEC.replace("units,1","units,NaN"), SPEC.replace("flatness,1e-8","flatness,0.1"),
            SPEC.replace("frame,0,0,0,1,0,0,0,1,0","frame,0,0,0,1,0,0,1,0,0"),
            SPEC.replace("bridge,69,0.5,0.5","bridge,69,1.1,0.5"),
            format!("{SPEC}fixed,1\n"), format!("{SPEC}bridge,69,0.2,0.2\n")] {
            assert!(import(OBJ,&spec).is_err());
        }
        assert!(import(&OBJ.replace("v 0.5 0.5 0","v 0.5 0.5 0.01"),SPEC).is_err());
        assert!(import(&format!("{OBJ}f 2 3 6\n"),SPEC).is_err());
    }
    #[test]
    fn region_sections_change_real_mass_without_shader_guessing() {
        let obj=OBJ.replace("f 4 5 6","usemtl maple\nf 4 5 6");
        let spec=format!("{SPEC}material,maple,0.008,750,1.2e10,1.2e9,0.3,9e8,0.2\n");
        let a=import(&obj,&spec).unwrap();
        let prepared=BoardGeometry::read(&a.fsb).unwrap().prepare(&[69],1000.0).unwrap();
        assert!((prepared.mass_kg-4.8).abs()<1e-10);
        assert!(import(&obj,SPEC).is_err());
    }
    #[test]
    fn native_steinway_sections_ribs_and_all_bridge_stations_round_trip() {
        let native=super::super::steinway_d::build(4).unwrap().geometry;
        let (obj,spec)=export(&native).unwrap();
        let imported=import(&obj,&spec).unwrap();
        let count=|s:&str,tag:&str|s.lines().filter(|line|line.starts_with(tag)).count();
        for tag in ["node,","triangle,","fixed,","stiffener,","bridge,"] {
            assert_eq!(count(&native,tag),count(&imported.fsb,tag),"{tag}");
        }
        assert_eq!(count(&imported.fsb,"bridge,"),88);
        for tag in ["triangle,","stiffener,"] {
            let parsed=|s:&str|s.lines().filter(|l|l.starts_with(tag))
                .map(|l|l.split(',').skip(1).map(|f|f.parse::<f64>().unwrap()).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            assert_eq!(parsed(&native),parsed(&imported.fsb),"{tag}");
        }
        let fixed = |s:&str| s.lines().filter_map(|r|r.strip_prefix("fixed,"))
            .map(|n|n.parse::<usize>().unwrap()).collect::<BTreeSet<_>>();
        assert_eq!(fixed(&native),fixed(&imported.fsb));
        BoardGeometry::read(&imported.fsb).unwrap();
    }
    #[test]
    fn beam_paths_are_remapped_and_missing_nodes_refuse() {
        let beam="stiffener,9e9,6e8,0.0005,1e-8,1e-8,-0.015,400,2,6,4\n";
        let imported=import(OBJ,&format!("{SPEC}{beam}")).unwrap();
        let row=imported.fsb.lines().find(|r|r.starts_with("stiffener,")).unwrap();
        assert!(row.ends_with(",0,4,2"));
        assert!(import(OBJ,&format!("{SPEC}{}",beam.replace(",2,6,4",",1,6,4"))).is_err());
        assert!(import(OBJ,&format!("{SPEC}{beam}{beam}")).is_err());
    }
}
