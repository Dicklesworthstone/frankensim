//! Actual planar FE boundary patches feeding the original compact-neck owner.
//! A reference Neumann face is a pressure-basis integration surface, not a
//! material closing wall. Its motion is an explicit uniform fluid port.
use super::*;
use crate::render::plate::impact::cavity::neck::CavityNeck;
use std::collections::BTreeSet;

/// Admission bound shared with the original compact-neck owner.
pub const MAX_SHELL_CAVITY_OPENINGS: usize = 8;

/// One fixed compact mouth on the original interior volume. Its pressure trace
/// and area are integrated from these faces, never supplied as modal numbers.
/// The internal length includes only the duct and any explicitly declared
/// interior correction. Resistance contains only retained internal loss. The
/// exterior model supplies its entire radiation impedance once.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellCavityOpening {
    pub boundary_faces: Vec<[u32; 3]>,
    pub internal_length_m: f64,
    pub internal_resistance_pa_s_m3: f64,
}

#[derive(Clone, Debug)]
pub(super) struct PreparedOpening {
    pub(super) source: ShellCavityOpening,
    /// Original coordinates in OUTWARD order, shared with the exterior owner.
    pub(super) triangles: Vec<[[f64; 3]; 3]>,
    area_m2: f64,
    pub(super) extent_m: f64,
}
impl PreparedOpening {
    pub(super) fn neck(&self, fe: &HelmholtzFe<'_>, shapes: &[Vec<f64>],
        gate: &CancelGate) -> Result<CavityNeck, ShellCavityError>
    {
        check_cancel(gate)?;
        let (load, area) = fe.boundary_load(&self.source.boundary_faces)
            .map_err(ShellCavityError::Acoustic)?;
        if (area-self.area_m2).abs() > 256.0*f64::EPSILON*self.area_m2 {
            return Err(bad("mouth quadrature and original physical area disagree"));
        }
        let mut averages = Vec::with_capacity(shapes.len());
        for shape in shapes {
            check_cancel(gate)?;
            let mean = load.iter().zip(shape).map(|(a,b)| a*b).sum::<f64>()/area;
            if !mean.is_finite() { return Err(bad("mouth pressure projection overflow")); }
            averages.push(mean);
        }
        Ok(CavityNeck {
            area_m2: area, effective_length_m: self.source.internal_length_m,
            resistance_pa_s_m3: self.source.internal_resistance_pa_s_m3,
            pressure_shape_averages: averages, initial_volume_m3: 0.0,
            initial_flow_m3_s: 0.0,
        })
    }
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    while parents[node] != node {
        parents[node] = parents[parents[node]];
        node = parents[node];
    }
    node
}

pub(super) fn prepare(parts: &[ShellPart<'_>], fe: &HelmholtzFe<'_>,
    sources: Vec<ShellCavityOpening>, gate: &CancelGate)
    -> Result<Vec<PreparedOpening>, ShellCavityError>
{
    if sources.len() > MAX_SHELL_CAVITY_OPENINGS {
        return Err(bad("compact cavity supports at most eight explicit mouths"));
    }
    for &(shell,rows) in parts {
        if rows.iter().any(|row|row.len() != 6*shell.mesh().nodes.len()) {
            return Err(bad("mouth rim admission requires the complete original shell displacement fields"));
        }
    }
    let mut selected = BTreeSet::new();
    let solid: BTreeSet<_> = parts.iter().flat_map(|(shell,_)|
        shell.mesh().tris.iter().map(|triangle| key(triangle.map(|n| shell.mesh().nodes[n])))).collect();
    let mut result = Vec::with_capacity(sources.len());
    for source in sources {
        check_cancel(gate)?;
        if source.boundary_faces.is_empty() || source.boundary_faces.len() > fe.boundary_face_count()
            || !source.internal_length_m.is_finite() || source.internal_length_m <= 0.0
            || !source.internal_resistance_pa_s_m3.is_finite() || source.internal_resistance_pa_s_m3 < 0.0
        { return Err(bad("mouth requires explicit boundary faces and finite positive internal length/passive loss")); }
        let mut triangles = Vec::with_capacity(source.boundary_faces.len());
        let mut area_m2 = 0.0;
        let mut area_centroid = [0.0;3];
        let mut first: Option<([f64;3],[f64;3])> = None;
        let mut parents: Vec<_> = (0..source.boundary_faces.len()).collect();
        let mut edges = BTreeMap::<[u32;2],(usize,usize)>::new();
        for (index,&face) in source.boundary_faces.iter().enumerate() {
            check_cancel(gate)?;
            let (mut points, normal, area) = fe.boundary_geometry(face)
                .map_err(ShellCavityError::Acoustic)?;
            let k = key(points);
            if !selected.insert(k) || solid.contains(&k) {
                return Err(bad("mouth repeats a physical face or covers an existing solid triangle"));
            }
            let (origin, reference) = *first.get_or_insert((points[0],normal));
            let alignment: f64 = normal.iter().zip(reference).map(|(a,b)|a*b).sum();
            let scale = points.iter().flatten().chain(origin.iter()).map(|v|v.abs())
                .fold(source.internal_length_m,f64::max);
            if alignment < 1.0-512.0*f64::EPSILON || points.iter().any(|p|
                (0..3).map(|i|(p[i]-origin[i])*reference[i]).sum::<f64>().abs()
                    > 512.0*f64::EPSILON*scale)
            { return Err(bad("one compact mouth must be an outward planar patch")); }
            for &(shell,rows) in parts {
                for (node,p) in shell.mesh().nodes.iter().enumerate() {
                    if points.contains(p) && rows.iter().any(|row| (0..3).any(|axis|row[6*node+axis] != 0.0)) {
                        return Err(bad("compact mouth rim must be fixed in every retained structural mode"));
                    }
                }
            }
            let a = [points[1][0]-points[0][0],points[1][1]-points[0][1],points[1][2]-points[0][2]];
            let b = [points[2][0]-points[0][0],points[2][1]-points[0][1],points[2][2]-points[0][2]];
            let cross = [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
            if cross.iter().zip(normal).map(|(a,b)|a*b).sum::<f64>() < 0.0 { points.swap(1,2); }
            triangles.push(points);
            area_m2 += area;
            for axis in 0..3 {
                area_centroid[axis] += area*(points[0][axis]+points[1][axis]+points[2][axis])/3.0;
            }
            for mut edge in [[face[0],face[1]],[face[1],face[2]],[face[2],face[0]]] {
                edge.sort_unstable();
                if let Some((other,count)) = edges.get_mut(&edge) {
                    *count += 1;
                    if *count > 2 { return Err(bad("mouth patch has a nonmanifold edge")); }
                    let a = root(&mut parents,index); let b = root(&mut parents,*other);
                    parents[a] = b;
                } else { edges.insert(edge,(index,1)); }
            }
        }
        let first = root(&mut parents,0);
        if (1..parents.len()).any(|i|root(&mut parents,i) != first) {
            return Err(bad("one compact mouth must be connected through complete edges"));
        }
        let mut rim = BTreeMap::<u32,Vec<u32>>::new();
        for (edge,(_,count)) in edges {
            if count == 1 {
                rim.entry(edge[0]).or_default().push(edge[1]);
                rim.entry(edge[1]).or_default().push(edge[0]);
            }
        }
        if rim.len() < 3 || rim.values().any(|neighbors|neighbors.len() != 2) {
            return Err(bad("compact mouth requires one simple boundary loop"));
        }
        let mut todo = vec![*rim.keys().next().expect("admitted mouth rim")];
        let mut seen = BTreeSet::new();
        while let Some(node) = todo.pop() {
            if seen.insert(node) { todo.extend(&rim[&node]); }
        }
        if seen.len() != rim.len() || !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(bad("compact mouth has holes or unrepresentable physical area"));
        }
        let centroid = area_centroid.map(|x|x/area_m2);
        let mut extent_m = source.internal_length_m;
        for point in triangles.iter().flatten() {
            extent_m = extent_m.max(det::hypot(
                det::hypot(point[0]-centroid[0],point[1]-centroid[1]),point[2]-centroid[2]));
        }
        if !extent_m.is_finite() { return Err(bad("mouth geometric extent overflow")); }
        result.push(PreparedOpening {source,triangles,area_m2,extent_m});
    }
    Ok(result)
}
