//! Cold import of actual rigid cabinet/lid parts, with explicit SI placement.
//! The moving board keeps its native source-facet map. Selected OBJ parts get
//! zero prescribed velocity and participate in that SAME exterior BEM solve.
//! This is static geometry: no hinges in the audio loop, collision repair,
//! elastic cabinet law, visual-material inference or new acoustic solver.
use super::{Boundary, SpherePanels, MAX_OBJ_BYTES, MAX_PANELS, MAX_SPEC_BYTES,
    closed_components, cross, dot, norm, number, sub};
use fs_math::det;
use std::{collections::BTreeMap, fmt::Write as _, io::Read, path::Path};

pub const HEADER: &str = "frankensim-piano-rigid-assembly-v1";
const MAX_PARTS: usize = 64;
type Triangle = [[f64; 3]; 3];

#[derive(Clone, Debug)]
struct Pose {
    pivot: [f64; 3], axis: [f64; 3], degrees: f64, translation: [f64; 3],
}
impl Pose {
    fn point(&self, raw: [f64; 3], scale: f64, origin: [f64; 3]) -> Result<[f64; 3], String> {
        let p = std::array::from_fn(|c| scale * (raw[c] - origin[c]) - self.pivot[c]);
        // A proper rotation preserves outward winding. Reflection and negative
        // unit scales are not admitted as a way to repair inward source meshes.
        let theta = self.degrees * std::f64::consts::PI / 180.;
        let (s, c) = (det::sin(theta), det::cos(theta));
        let skew = cross(self.axis, p); let along = dot(self.axis, p);
        let out: [f64; 3] = std::array::from_fn(|i| self.translation[i] + self.pivot[i]
            + c * p[i] + s * skew[i] + (1. - c) * along * self.axis[i]);
        if out.iter().any(|x| !x.is_finite() || x.abs() > 100.) {
            return Err("posed rigid geometry exceeds the finite 100 metre scene budget".into());
        }
        Ok(out)
    }
}
#[derive(Clone, Debug)]
struct Selection { path: String, label: String, scale: f64, origin: [f64; 3] }
#[derive(Debug)]
struct Part {
    name: String, selection: Selection, pose: Pose, triangles: Vec<Triangle>, excluded: usize,
}

/// Immutable admitted coordinates, not file references reopened during a render.
/// MTLs, textures, network URLs and other nested asset references are never read.
#[derive(Debug)]
pub struct Assembly { source: String, parts: Vec<Part>, panels: usize }

fn read(path: &Path, cap: usize) -> Result<String, String> {
    let mut bytes = Vec::new();
    std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?
        .take(cap as u64 + 1).read_to_end(&mut bytes).map_err(|e| e.to_string())?;
    if bytes.len() > cap { return Err(format!("{} exceeds the rigid-asset byte budget", path.display())); }
    String::from_utf8(bytes).map_err(|_| format!("{} must be UTF-8 OBJ text", path.display()))
}
fn alias(name: &str) -> bool {
    !name.is_empty() && name.len() <= 128 && name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
}
fn vector(f: &[&str]) -> Result<[f64; 3], String> { Ok([number(f[0])?, number(f[1])?, number(f[2])?]) }

impl Assembly {
    /// Relative OBJ paths resolve against the assembly file's directory, not
    /// the process working directory. Only explicitly listed OBJ files open.
    pub fn load(path: &str) -> Result<Self, String> {
        let path = Path::new(path); let text = read(path, MAX_SPEC_BYTES)?;
        let base = path.parent().unwrap_or_else(|| Path::new("."));
        Self::from_text(&text, |name| read(&base.join(name), MAX_OBJ_BYTES))
    }
    /// Every part needs one source selection and one pose. The caller owns
    /// source acquisition; this hook also allows assets already held in memory.
    /// Total unique source text is capped at 32 MiB; complete selected geometry
    /// must fit the SAME 2048-panel limit as the acoustic scene.
    pub fn from_text(text: &str, mut load: impl FnMut(&str) -> Result<String, String>) -> Result<Self, String> {
        if text.len() > MAX_SPEC_BYTES { return Err("rigid assembly exceeds 64 KiB".into()); }
        let mut rows = text.lines().map(str::trim).filter(|r| !r.is_empty() && !r.starts_with('#'));
        if rows.next() != Some(HEADER) { return Err(format!("expected {HEADER}")); }
        let mut source = None; let mut selections = BTreeMap::new(); let mut poses = BTreeMap::new();
        for row in rows {
            let f: Vec<_> = row.split(',').map(str::trim).collect();
            match f[0] {
                "source" if f.len() >= 3 => {
                    if source.is_some() || !["estimated", "mixed", "published", "measured"].contains(&f[1])
                        || f[2..].iter().all(|x| x.is_empty()) {
                        return Err("rigid assembly requires one attributed source row".into());
                    }
                    source = Some(f[1..].join(","));
                }
                "part" if f.len() == 8 => {
                    if !alias(f[1]) || f[2].is_empty() || f[3].is_empty() || selections.len() == MAX_PARTS {
                        return Err("part needs an alias, explicit OBJ path and exact object/group label".into());
                    }
                    let scale = number(f[4])?;
                    if !(1e-9..=1e6).contains(&scale) { return Err("rigid OBJ unit scale must be positive".into()); }
                    let selection = Selection { path: f[2].into(), label: f[3].into(), scale, origin: vector(&f[5..8])? };
                    if selections.insert(f[1].to_owned(), selection).is_some() { return Err("duplicate rigid part alias".into()); }
                }
                "pose" if f.len() == 12 => {
                    let pose = Pose { pivot: vector(&f[2..5])?, axis: vector(&f[5..8])?,
                        degrees: number(f[8])?, translation: vector(&f[9..12])? };
                    if !alias(f[1]) || poses.len() == MAX_PARTS || (norm(pose.axis) - 1.).abs() > 1e-12
                        || pose.degrees.abs() > 360. || pose.pivot.iter().chain(&pose.translation).any(|x| x.abs() > 100.) {
                        return Err("pose needs an SI pivot/translation, unit axis and signed angle within 360 degrees".into());
                    }
                    if poses.insert(f[1].to_owned(), pose).is_some() { return Err("duplicate rigid pose".into()); }
                }
                _ => return Err(format!("unknown rigid assembly row or field count: {}", f[0])),
            }
        }
        if source.is_none() || selections.is_empty() || selections.len() != poses.len()
            || selections.keys().any(|k| !poses.contains_key(k)) {
            return Err("every rigid selection requires exactly one pose and source attribution".into());
        }
        let mut documents = BTreeMap::new(); let mut bytes = 0usize; let mut panels = 0usize; let mut parts = Vec::new();
        for (name, selection) in selections {
            if !documents.contains_key(&selection.path) {
                let text = load(&selection.path)?;
                bytes = bytes.checked_add(text.len()).filter(|n| *n <= MAX_OBJ_BYTES)
                    .ok_or("unique rigid OBJ inputs exceed 32 MiB")?;
                let document = fs_io::obj::read_obj_document(&text).map_err(|e| format!("{}: {e}", selection.path))?;
                documents.insert(selection.path.clone(), document);
            }
            let doc = &documents[&selection.path];
            let pose = poses.remove(&name).ok_or("missing admitted rigid pose")?;
            let mut triangles = Vec::new();
            // Transform each selected source vertex once, preserving exact
            // shared coordinates and seams even when several groups meet.
            let mut points = BTreeMap::new();
            for region in &doc.regions {
                if !region.has_label(&selection.label) { continue; }
                for t in &doc.soup.triangles[region.triangles.clone()] {
                    if panels == MAX_PANELS { return Err("selected rigid parts exceed the 2048-panel scene budget".into()); }
                    let mut triangle = [[0.; 3]; 3];
                    for (i, &index) in t.iter().enumerate() {
                        let p = if let Some(&p) = points.get(&index) { p } else {
                            let p = doc.soup.positions[index as usize];
                            let p = pose.point([p.x, p.y, p.z], selection.scale, selection.origin)?;
                            points.insert(index, p); p
                        };
                        triangle[i] = p;
                    }
                    triangles.push(triangle); panels += 1;
                }
            }
            if triangles.is_empty() { return Err(format!("{name}: no faces match exact label {:?}", selection.label)); }
            // Reuse the existing triangle and closed-component admission; never
            // cap an open artwork mesh or flip a component to force acceptance.
            SpherePanels::from_triangles(triangles.clone()).map_err(|e| format!("{name}: {e}"))?;
            closed_components(&triangles).map_err(|e| format!("{name}: {e}"))?;
            let excluded = doc.soup.triangles.len() - triangles.len();
            parts.push(Part { name, selection, pose, triangles, excluded });
        }
        let combined: Vec<_> = parts.iter().flat_map(|p| p.triangles.iter().copied()).collect();
        closed_components(&combined)?;
        Ok(Self { source: source.ok_or("missing source")?, parts, panels })
    }
    pub fn panel_count(&self) -> usize { self.panels }
    pub fn report(&self) -> String {
        let mut out = format!("Rigid assembly [{}]: {} parts / {} panels; static acoustically rigid scatterers, not structural masses or damping.",
            self.source, self.parts.len(), self.panels);
        for p in &self.parts {
            let _ = write!(out, " {}: {:?} label {:?}, {} selected / {} excluded, scale {} m/unit, origin {:?}, pivot {:?}, axis {:?}, {} deg, translation {:?} m;",
                p.name, p.selection.path, p.selection.label, p.triangles.len(), p.excluded, p.selection.scale,
                p.selection.origin, p.pose.pivot, p.pose.axis, p.pose.degrees, p.pose.translation);
        }
        out
    }
    /// Append without touching the board's known embeddings, modes or motion.
    /// Zero added velocities impose the rigid boundary in the joint BEM; these
    /// panels are not extra sources, an IR, or a post-render cabinet effect.
    pub fn attach(&self, mut body: Boundary) -> Result<Boundary, String> {
        let original = body.surface.triangles().ok_or("rigid composition needs retained source triangles")?;
        let count = original.len().checked_add(self.panels).filter(|n| *n <= MAX_PANELS)
            .ok_or("combined board and rigid parts exceed the 2048-panel BEM budget")?;
        if body.weights.is_empty() || body.weights.iter().any(|r| r.len() != original.len() || r.iter().any(|v| !v.is_finite())) {
            return Err("rigid composition needs every finite board velocity row".into());
        }
        let mut triangles = Vec::with_capacity(count); triangles.extend_from_slice(original);
        for part in &self.parts { triangles.extend_from_slice(&part.triangles); }
        let components = closed_components(&triangles)?;
        let center = std::array::from_fn(|c| {
            let lo = triangles.iter().flatten().map(|p| p[c]).fold(f64::INFINITY, f64::min);
            let hi = triangles.iter().flatten().map(|p| p[c]).fold(f64::NEG_INFINITY, f64::max);
            0.5 * (lo + hi)
        });
        let radius = triangles.iter().flatten().map(|&p| norm(sub(p, center))).fold(0.0_f64, f64::max);
        let surface = SpherePanels::from_triangles(triangles).map_err(|e| e.to_string())?;
        for row in &mut body.weights { row.resize(count, 0.); }
        Ok(Boundary { surface, weights: body.weights, center, radius, components })
    }
    /// Deterministic posed SI geometry for inspection; negative face indices
    /// make this safe to append to another valid OBJ. No material files copied.
    pub fn obj(&self) -> String {
        let mut text = format!("# {}\n# Geometry only; all exported parts are acoustically rigid.\n", self.report());
        for part in &self.parts {
            let _ = writeln!(text, "o {}", part.name);
            for triangle in &part.triangles {
                for p in triangle { let _ = writeln!(text, "v {:.17e} {:.17e} {:.17e}", p[0], p[1], p[2]); }
                text.push_str("f -3 -2 -1\n");
            }
        }
        text
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{tests as fixtures, Specification, MotionSurface, TAU};
    fn input(angle: f64) -> String {
        format!("{HEADER}\nsource,estimated,authored unit/hinge regression\npart,lid,parts.obj,Lid wood,0.001,0,0,0\npose,lid,0,0,0.06,1,0,0,{angle},0,0,0\n")
    }
    fn asset() -> String { fixtures::box_obj("Lid wood", [0.,0.,60.], [100.,100.,20.]) }
    fn motion() -> MotionSurface {
        MotionSurface::new(fs_plate::ShellMesh::new(vec![[0.,0.,0.],[0.1,0.,0.],[0.1,0.1,0.],[0.,0.1,0.]],
            vec![[0,1,2],[0,2,3]]).unwrap(), vec![vec![[0.,0.,1.,0.,0.,0.];4]]).unwrap()
    }
    #[test]
    fn selected_parts_convert_units_and_rotate_about_the_actual_hinge_without_reading_materials() {
        let text = format!("{}part,copy,parts.obj,Lid wood,0.001,0,0,0\npose,copy,0,0,0.06,1,0,0,0,1,0,0\n", input(90.));
        let mut reads = 0;
        let assembled = Assembly::from_text(&text, |path| { reads += 1; assert_eq!(path,"parts.obj");
            Ok(format!("mtllib not-opened.mtl\n{}{}",asset(),fixtures::box_obj("Excluded",[0.;3],[1.;3]))) }).unwrap();
        assert_eq!(reads,1); assert_eq!(assembled.panel_count(),24);
        let lid = assembled.parts.iter().find(|p|p.name=="lid").unwrap();
        assert_eq!(lid.excluded,12);
        for p in lid.triangles.iter().flatten() {
            assert!(p[1]>=-0.0200000000001 && p[1]<=1e-12);
            assert!(p[2]>=0.0599999999999 && p[2]<=0.160000000001);
        }
        let obj = assembled.obj();
        let parsed = fs_io::obj::read_obj_document(&obj).unwrap();
        assert_eq!(parsed.soup.triangles.len(),24); assert!(obj.contains("o lid\n"));
        assert_eq!(obj,assembled.obj());
    }
    #[test]
    fn posed_lid_changes_the_actual_pressure_but_not_the_board_motion_map() {
        let spec = Specification::read(&fixtures::specification()).unwrap();
        let board = || Boundary::from_obj(&fixtures::box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]),&spec,&motion()).unwrap();
        let a = Assembly::from_text(&input(0.), |_| Ok(asset())).unwrap().attach(board()).unwrap();
        let b = Assembly::from_text(&input(45.), |_| Ok(asset())).unwrap().attach(board()).unwrap();
        let bare = board();
        for posed in [&a,&b] {
            assert_eq!(posed.components,2); assert_eq!(&posed.weights[0][..12],&bare.weights[0]);
            assert!(posed.weights[0][12..].iter().all(|x|*x==0.));
            assert_eq!(&posed.surface.triangles().unwrap()[..12],bare.surface.triangles().unwrap());
        }
        let w=[TAU*100.,TAU*200.];
        let pa=a.sample_grid(&w,&spec.receivers,spec.medium,6.).unwrap();
        let pb=b.sample_grid(&w,&spec.receivers,spec.medium,6.).unwrap();
        assert!((pa.values[0][0][1]-pb.values[0][0][1]).abs()>1e-6*pa.values[0][0][1].abs());
        assert_ne!(a.center,b.center); // receiver flight guard must use the posed whole body
    }
    #[test]
    fn incomplete_pose_and_unusable_source_geometry_do_not_create_an_assembly() {
        for bad in [input(0.).replace("pose,lid,0,0,0.06,1,0,0,0,0,0,0\n", ""),
            input(0.).replace("1,0,0,0,0,0,0", "2,0,0,0,0,0,0"),
            input(0.).replace("0.001", "-0.001"), input(0.).replace("0.001", "NaN"),
            format!("{}pose,unknown,0,0,0,1,0,0,0,0,0,0\n",input(0.))] {
            assert!(Assembly::from_text(&bad, |_| panic!("malformed manifest must precede asset I/O")).is_err());
        }
        assert!(Assembly::from_text(&input(0.).replace("Lid wood","lid wood"), |_| Ok(asset())).is_err());
        let open=asset().replace("f -8 -6 -7\n","");
        assert!(Assembly::from_text(&input(0.), |_| Ok(open.clone())).is_err());
        let a=Assembly::from_text(&input(0.), |_| Ok(asset())).unwrap();
        let too_many=SpherePanels::from_triangles(vec![a.parts[0].triangles[0];MAX_PANELS]).unwrap();
        assert!(a.attach(Boundary {surface:too_many,weights:vec![vec![0.;MAX_PANELS]],center:[0.;3],radius:1.,components:1}).is_err());
    }
}
