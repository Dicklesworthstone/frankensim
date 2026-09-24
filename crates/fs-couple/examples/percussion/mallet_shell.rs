//! A flat felt face above the actual curved, finite-thickness upper skin.
//! Contact gaps come from the complete circular footprint, not its centroid.
//! Four positive-area sites retain separate strain histories and rotational work.
use super::{Error, Spec, Stroke, PadSide, PadSite, FeltStriker};
use super::super::playing;
use fs_plate::shell::{survey::MeshShell, reduction::{ShellReduction, radiation::ShellFace}};
use std::collections::BTreeMap;

impl Spec {
    #[allow(clippy::too_many_arguments)]
    pub fn compile_shell(&self,reduction:&ShellReduction,shell:&MeshShell,stroke:Stroke,
        coordinate:usize,first_shell:usize,total:usize)->Result<FeltStriker,Error>
    {
        let end=first_shell.checked_add(reduction.mode_count()).ok_or("mallet shell address overflow")?;
        if total>fs_couple::render::plate::impact::MAX_IMPACT_MODES || end>total || coordinate>=total
            || (first_shell..end).contains(&coordinate) || !stroke.speed_m_s.is_finite()
            || !(0.0..=20.0).contains(&stroke.speed_m_s) {
            return Err("mallet needs a separate striker address and complete bounded shell basis".into());
        }
        let center=stroke.position_m.ok_or("cymbal mallets require an explicit strike XY station")?;
        let nodes=reduction.surface_positions(&shell.nodal_thickness_m,ShellFace::Positive)?;
        let mut plane=contact_plane(&nodes,&shell.mesh.tris,center,self.radius_m)?;
        let mut ports=Vec::with_capacity(4);
        for point in self.points(center) {
            // Locate on the actual positive SKIN chart, not the midsurface or
            // the hi-hat's inner face. Original triangle order is retained.
            let (face,bary)=playing::shell_location(&nodes,&shell.mesh.tris,point)?;
            let port=reduction.surface_point_port(&shell.nodal_thickness_m,face,bary,
                ShellFace::Positive,[0.,0.,1.])?;
            plane=plane.max(port.position_m[2]); // include sampled roundoff, never a negative gap
            ports.push(port);
        }
        let area=std::f64::consts::PI*self.radius_m*self.radius_m/4.0;
        let mut sites=Vec::with_capacity(4);let mut clearances=Vec::with_capacity(4);
        for port in ports {
            let mut weights=vec![0.;total];weights[first_shell..end].copy_from_slice(&port.weights);
            sites.push(PadSite {weights,area_m2:area});clearances.push(plane-port.position_m[2]);
        }
        let mut jaw=self.jaw.clone();jaw.side=PadSide::Positive;jaw.initial_velocity_m_s=stroke.speed_m_s;
        // Positive shell travel is UP, positive mallet travel/force is DOWN.
        // Area belongs to the horizontal felt face, not the sloped shell facet.
        let tip=FeltStriker::new_with_clearances(total,coordinate,&sites,&jaw,&clearances)?;
        eprintln!("curved felt strike: radius_m={}, contact_plane_z_m={plane}, max_extra_gap_m={}, sites=4; flat axial face, independent strain histories, not a calibrated mallet",
            self.radius_m,clearances.iter().copied().fold(0.0_f64,f64::max));
        Ok(tip)
    }
}
fn cross(a:[f64;2],b:[f64;2])->f64{a[0]*b[1]-a[1]*b[0]}
fn sub(a:[f64;2],b:[f64;2])->[f64;2]{[a[0]-b[0],a[1]-b[1]]}
fn dot(a:[f64;2],b:[f64;2])->f64{a[0]*b[0]+a[1]*b[1]}
fn distance(a:[f64;2],b:[f64;2],p:[f64;2])->f64 {
    let d=sub(b,a);let q=sub(p,a);let t=(dot(q,d)/dot(d,d)).clamp(0.,1.);
    (q[0]-t*d[0]).hypot(q[1]-t*d[1])
}

/// P1-height maximum over the WHOLE disk. The chart must face upward everywhere,
/// have unit winding over this disk, and keep every outer/hole edge outside it.
/// These checks reject clipping, holes, projected folds and multiple coverage.
/// They are floating-point geometry admission, not exact-predicate certification.
fn contact_plane(nodes:&[[f64;3]],triangles:&[[usize;3]],center:[f64;2],radius:f64)->Result<f64,Error> {
    if nodes.is_empty() || triangles.is_empty() || center.iter().any(|x|!x.is_finite())
        || !radius.is_finite() || radius<=0. || nodes.iter().flatten().any(|x|!x.is_finite()) {
        return Err("mallet footprint requires finite nonempty skin geometry and positive radius".into());
    }
    let mut edges=BTreeMap::<(usize,usize),(usize,usize,usize)>::new();
    let mut peak=f64::NEG_INFINITY;
    for tri in triangles {
        if tri.iter().any(|&i|i>=nodes.len()) {return Err("mallet skin connectivity mismatch".into());}
        let t=tri.map(|i|nodes[i]);
        if let Some(h)=triangle_peak(t,center,radius)? {peak=peak.max(h);}
        for j in 0..3 {
            let (a,b)=(tri[j],tri[(j+1)%3]);
            let e=edges.entry((a.min(b),a.max(b))).or_insert((a,b,0));
            if e.2>1 || (e.2==1 && (e.0,e.1)!=(b,a)) {
                return Err("mallet skin must have consistently oriented manifold edges".into());
            }
            e.2+=1;
        }
    }
    let mut winding=0_i64;
    for &(a,b,count) in edges.values() {
        if count!=1 {continue;}
        let a=[nodes[a][0],nodes[a][1]];let b=[nodes[b][0],nodes[b][1]];
        let d=distance(a,b,center);
        if !d.is_finite() || d<=radius {
            return Err("full mallet disk crosses a shell rim or mounting hole; no clipping or snapping".into());
        }
        let side=cross(sub(b,a),sub(center,a));
        if !side.is_finite() {return Err("mallet boundary winding overflows".into());}
        if a[1]<=center[1] && b[1]>center[1] && side>0. {winding+=1;}
        if b[1]<=center[1] && a[1]>center[1] && side<0. {winding-=1;}
    }
    if winding!=1 || !peak.is_finite() {return Err("mallet footprint must cover exactly one upward-facing skin sheet".into());}
    playing::shell_location(nodes,triangles,center)?;
    Ok(peak)
}

// Maximize an affine height on triangle intersect disk. Candidates are triangle
// vertices in the disk, edge/circle intersections, and the unconstrained circle
// support point in the height-gradient direction. No coarse perimeter sampling.
fn triangle_peak(t:[[f64;3];3],center:[f64;2],radius:f64)->Result<Option<f64>,Error> {
    let p=t.map(|p|[(p[0]-center[0])/radius,(p[1]-center[1])/radius]);
    let u=sub(p[1],p[0]);let v=sub(p[2],p[0]);let det=cross(u,v);
    if p.iter().flatten().any(|x|!x.is_finite()) || !det.is_finite() || det<=0. {
        return Err("flat cymbal mallet requires an upward, nonfolded XY skin chart".into());
    }
    let dz=[t[1][2]-t[0][2],t[2][2]-t[0][2]];
    let g=[(dz[0]*v[1]-dz[1]*u[1])/det,(u[0]*dz[1]-v[0]*dz[0])/det];
    if g.iter().any(|x|!x.is_finite()) {return Err("mallet skin height gradient overflows".into());}
    let inside=|q:[f64;2]| {
        let d=sub(q,p[0]);let a=cross(d,v)/det;let b=cross(u,d)/det;
        a>=-64.*f64::EPSILON && b>=-64.*f64::EPSILON && a+b<=1.+64.*f64::EPSILON
    };
    let mut peak=f64::NEG_INFINITY;
    let mut candidate=|q:[f64;2]| {peak=peak.max(t[0][2]+dot(g,sub(q,p[0])));};
    for &q in &p {if dot(q,q)<=1. {candidate(q);}}
    if inside([0.;2]) {candidate([0.;2]);}
    let norm=g[0].hypot(g[1]);
    if norm>0. {
        let q=[g[0]/norm,g[1]/norm];if inside(q){candidate(q);}
    }
    for j in 0..3 {
        let a=p[j];let d=sub(p[(j+1)%3],a);let dd=dot(d,d);
        if !dd.is_finite() || dd<=0. {return Err("mallet projected edge is unresolved".into());}
        let middle=-dot(a,d)/dd;
        let near=[a[0]+middle*d[0],a[1]+middle*d[1]];
        let square=1.-dot(near,near);
        if square<0. {continue;}
        let delta=(square/dd).sqrt();
        for r in [middle-delta,middle+delta] {
            if (0.0..=1.0).contains(&r) {candidate([a[0]+r*d[0],a[1]+r*d[1]]);}
        }
    }
    if peak==f64::NEG_INFINITY {return Ok(None);}
    if !peak.is_finite(){return Err("mallet reference contact plane overflows".into());}
    Ok(Some(peak))
}

#[cfg(test)]
#[path="mallet_shell_tests.rs"]
mod tests;
