//! Conservative finite-aperture Neumann data on the existing exterior BEM mesh.
//! The neck's OWN volume-flow coordinate supplies flux, never interior pressure.
//! Local conforming subdivision resolves the mouth without remeshing the heads.
//! Geometry remains the original piecewise-planar cylinder: refine its parent
//! azimuths separately. This is a prescribed-flow, one-way radiation image.
use super::*;
use fs_couple::render::plate::impact::cavity::{cylinder::SidewallAperture, neck::NeckRadiationPort};
use std::collections::BTreeSet;

#[derive(Clone)]
struct Face {
    nodes: [usize; 3],
    // Unrolled sidewall coordinates relative to this opening. None on heads/rims.
    chart: Option<[[f64; 2]; 3]>,
    parent: usize,
}
fn wrap(angle: f64) -> f64 {
    (angle + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
}
fn side_chart(t: &[[f64; 3]; 3], radius: f64, depth: f64, a: SidewallAperture)
    -> Option<[[f64; 2]; 3]>
{
    if t.iter().any(|p| (p[0].hypot(p[1])-radius).abs() > 1e-10*radius) { return None; }
    let origin = a.azimuth_rad.rem_euclid(std::f64::consts::TAU);
    let first = wrap(t[0][1].atan2(t[0][0])-origin);
    let mut theta = t.map(|p| first+wrap(p[1].atan2(p[0])-origin-first));
    // A triangle across the opposite azimuthal seam must not span the mouth.
    let center = theta.iter().sum::<f64>()/3.0;
    let shift = center-wrap(center);
    for value in &mut theta { *value -= shift; }
    Some(core::array::from_fn(|i| [radius*theta[i], depth/2.0-t[i][2]-a.axial_position_m]))
}
fn area(t: [[f64; 3]; 3]) -> f64 {
    let u: [f64; 3] = core::array::from_fn(|i| t[1][i]-t[0][i]);
    let v: [f64; 3] = core::array::from_fn(|i| t[2][i]-t[0][i]);
    let c = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]];
    0.5*c[0].hypot(c[1]).hypot(c[2])
}
fn near_mouth(t: [[f64; 2]; 3], radius: f64) -> bool {
    // Conservative box/circle intersection; surplus refinement is harmless.
    let distance: [f64; 2] = core::array::from_fn(|axis| {
        let low = t.iter().map(|p| p[axis]).fold(f64::INFINITY, f64::min);
        let high = t.iter().map(|p| p[axis]).fold(f64::NEG_INFINITY, f64::max);
        if low > 0.0 { low } else if high < 0.0 { -high } else { 0.0 }
    });
    distance[0].hypot(distance[1]) <= radius
}
fn longest(t: [[f64; 2]; 3]) -> f64 {
    (0..3).map(|i| (t[i][0]-t[(i+1)%3][0]).hypot(t[i][1]-t[(i+1)%3][1]))
        .fold(0.0, f64::max)
}
fn edge(a: usize, b: usize) -> (usize, usize) { (a.min(b), a.max(b)) }
fn contains(t: [[f64; 2]; 3], p: [f64; 2]) -> bool {
    let cross = |a: [f64; 2], b: [f64; 2], c: [f64; 2]|
        (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
    let determinant = cross(t[0],t[1],t[2]);
    if !determinant.is_finite() || determinant == 0.0 { return false; }
    let b = [cross(t[1],t[2],p)/determinant, cross(t[2],t[0],p)/determinant,
        cross(t[0],t[1],p)/determinant];
    b.iter().all(|v| v.is_finite() && *v >= -64.0*f64::EPSILON && *v <= 1.0+64.0*f64::EPSILON)
}

/// Split every marked edge on BOTH incident triangles (red/green refinement).
/// No hanging nodes, duplicated panels, changed parent velocities or smoothing.
fn refine(points: &mut Vec<[f64; 3]>, faces: &[Face], marked: &BTreeSet<(usize, usize)>,
    maximum: usize) -> Result<Vec<Face>, Error>
{
    let mut mids = BTreeMap::new();
    for &(a,b) in marked {
        mids.insert((a,b),points.len());
        points.push(core::array::from_fn(|i| points[a][i]/2.0+points[b][i]/2.0));
    }
    let mut out = Vec::new();
    for face in faces {
        let n = face.nodes;
        let m: [Option<usize>; 3] = core::array::from_fn(|i| mids.get(&edge(n[i],n[(i+1)%3])).copied());
        let count = m.iter().flatten().count();
        if out.len()+count+1 > maximum { return Err("aperture refinement exceeds the boundary panel budget".into()); }
        let ids = [n[0],n[1],n[2],m[0].unwrap_or(0),m[1].unwrap_or(0),m[2].unwrap_or(0)];
        let chart = face.chart.map(|t| [t[0],t[1],t[2],
            [(t[0][0]+t[1][0])/2.0,(t[0][1]+t[1][1])/2.0],
            [(t[1][0]+t[2][0])/2.0,(t[1][1]+t[2][1])/2.0],
            [(t[2][0]+t[0][0])/2.0,(t[2][1]+t[0][1])/2.0]]);
        let mut emit = |local: [usize; 3]| out.push(Face {
            nodes: local.map(|i| ids[i]), chart: chart.map(|c| local.map(|i| c[i])), parent: face.parent,
        });
        match count {
            0 => emit([0,1,2]),
            1 => {
                let a = m.iter().position(Option::is_some).unwrap(); let b=(a+1)%3; let c=(a+2)%3;
                emit([a,3+a,c]); emit([3+a,b,c]);
            }
            2 => {
                let missing = m.iter().position(Option::is_none).unwrap();
                let a=(missing+1)%3; let b=(a+1)%3; let c=(a+2)%3;
                emit([b,3+b,3+a]); emit([a,3+a,3+b]); emit([a,3+b,c]);
            }
            _ => { emit([0,3,5]); emit([3,1,4]); emit([5,4,2]); emit([3,4,5]); }
        }
    }
    Ok(out)
}

impl Boundary {
    /// Add the outward volume-flow port Q=b*p as a resolved sidewall source.
    /// The computational boundary stays closed: its mouth cap carries Neumann
    /// flux, not a rigid-wall condition. All original head velocities are kept.
    ///
    /// Positive-area quadrature matches the interior's unrolled circular chart.
    /// Each quadrature atom is assigned ONCE to its containing refined panel;
    /// panel velocity times physical area sums to b, including the curved-to-
    /// planar area conversion. No pressure-to-audio shortcut or extra 1/r gain.
    ///
    /// Up to 12 conforming refinements resolve chart edges to a/2 near the mouth.
    /// maximum_terms bounds quadrature-point/panel searches. max_panels cannot
    /// exceed the existing BEM ceiling. Neither bound is a convergence claim.
    pub fn with_unrolled_sidewall_aperture(self, radius: f64, depth: f64, opening: SidewallAperture,
        volume_port: NeckRadiationPort, max_panels: usize, gate: &CancelGate) -> Result<Self, Error>
    {
        if gate.is_requested() { return Err("aperture preparation cancelled".into()); }
        let a = opening.radius_m; let coordinate=volume_port.coordinate;
        let weight=volume_port.volume_weight_m2_per_sqrt_kg;
        let nominal_area=std::f64::consts::PI*a*a;
        if !radius.is_finite() || radius<=0.0 || !depth.is_finite() || depth<=0.0
            || !a.is_finite() || a<=0.0 || a>radius/10.0
            || !volume_port.area_m2.is_finite() || volume_port.area_m2<=0.0
            || (nominal_area-volume_port.area_m2).abs()>1e-12*volume_port.area_m2
            || !volume_port.effective_length_m.is_finite() || volume_port.effective_length_m<=0.0
            || std::f64::consts::TAU*1640.0/Medium::air().sound_speed*a.max(volume_port.effective_length_m)>0.3
            || !opening.azimuth_rad.is_finite() || !opening.axial_position_m.is_finite()
            || opening.axial_position_m-a<0.0 || opening.axial_position_m+a>depth
            || !(1..=64).contains(&opening.radial_rings) || !(8..=256).contains(&opening.angular_points)
            || !weight.is_finite() || weight<=0.0 || coordinate>=fs_couple::render::plate::impact::MAX_IMPACT_MODES
            || self.state_modes.contains(&coordinate) || self.state_modes.len()>=MAX_INPUTS
            || self.triangles.is_empty() || self.triangles.len()>max_panels || max_panels>MAX_PANELS
            || self.weights.len()!=self.state_modes.len()
            || self.weights.iter().any(|r| r.len()!=self.triangles.len() || r.iter().any(|v| !v.is_finite()))
            || self.triangles.iter().flatten().flatten().any(|v| !v.is_finite())
        { return Err("aperture requires finite interior sidewall geometry, an independent flow port and bounded data".into()); }
        let samples=opening.radial_rings*opening.angular_points;
        let mut points=Vec::new(); let mut lookup=BTreeMap::new(); let mut faces=Vec::new();
        for (parent,t) in self.triangles.iter().enumerate() {
            let nodes=t.map(|p| *lookup.entry(p.map(|x| if x==0.0 {0}else{x.to_bits()})).or_insert_with(|| {
                let index=points.len(); points.push(p); index
            }));
            if !area(*t).is_finite() || area(*t)<=0.0 { return Err("aperture boundary has a degenerate panel".into()); }
            faces.push(Face {nodes,chart:side_chart(t,radius,depth,opening),parent});
        }
        check_closed(&faces.iter().map(|f| f.nodes).collect::<Vec<_>>())?;
        let mut resolved=false;
        for level in 0..=12 {
            if gate.is_requested() { return Err("aperture preparation cancelled".into()); }
            let mut marked=BTreeSet::new();
            for f in &faces {
                if f.chart.is_some_and(|t| near_mouth(t,a) && longest(t)>a/2.0) {
                    for i in 0..3 { marked.insert(edge(f.nodes[i],f.nodes[(i+1)%3])); }
                }
            }
            if marked.is_empty() { resolved=true; break; }
            if level==12 { break; }
            if samples.checked_mul(faces.len()).is_none_or(|n| n>opening.maximum_terms) {
                return Err("aperture projection exceeds its quadrature work budget".into());
            }
            faces=refine(&mut points,&faces,&marked,max_panels)?;
        }
        if !resolved { return Err("aperture spatial refinement exhausted; refine the source cylinder".into()); }
        if samples.checked_mul(faces.len()).is_none_or(|n| n>opening.maximum_terms) {
            return Err("aperture projection exceeds its quadrature work budget".into());
        }
        check_closed(&faces.iter().map(|f| f.nodes).collect::<Vec<_>>())?;
        let mut counts=vec![0usize;faces.len()];
        for ring in 0..opening.radial_rings {
            let r=a*((ring as f64+0.5)/opening.radial_rings as f64).sqrt();
            for sample in 0..opening.angular_points {
                if gate.is_requested() { return Err("aperture projection cancelled".into()); }
                let phi=std::f64::consts::TAU*(sample as f64+0.5)/opening.angular_points as f64;
                let p=[r*det::cos(phi),r*det::sin(phi)];
                let index=faces.iter().position(|f| f.chart.is_some_and(|t| contains(t,p)))
                    .ok_or("aperture quadrature point has no containing sidewall panel")?;
                counts[index]+=1;
            }
        }
        let triangles: Vec<_>=faces.iter().map(|f| f.nodes.map(|i| points[i])).collect();
        let mut weights: Vec<Vec<f64>>=self.weights.iter()
            .map(|row| faces.iter().map(|f| row[f.parent]).collect()).collect();
        let mut flow=vec![0.0;faces.len()];
        for (i,&count) in counts.iter().enumerate() {
            let panel_area=area(triangles[i]);
            if !panel_area.is_finite() || panel_area<=0.0 { return Err("refined aperture panel is degenerate".into()); }
            flow[i]=weight*(count as f64/samples as f64)/panel_area;
            if !flow[i].is_finite() || (count>0 && flow[i]<=0.0) { return Err("aperture normal velocity is unrepresentable".into()); }
        }
        weights.push(flow); let mut state_modes=self.state_modes; state_modes.push(coordinate);
        if gate.is_requested() { return Err("aperture preparation cancelled".into()); }
        Ok(Self {triangles,weights,state_modes})
    }
}

#[cfg(test)]
#[path = "curved_aperture_tests.rs"]
mod tests;

/// Uniform cold refinement of the SAME polyhedral boundary and piecewise-constant
/// source fields. Reuse the conforming aperture refiner; never smooth/project a
/// midpoint onto a different surface or average unrelated source-mode rows.
pub(super) fn uniform_refinement(boundary: &Boundary, levels: u32, maximum: usize,
    gate: &CancelGate) -> Result<Boundary, Error>
{
    if gate.is_requested() { return Err("radiation refinement cancelled".into()); }
    let count=boundary.triangles.len();
    if levels>4 || count==0 || maximum>fs_bem::helmholtz::MAX_DENSE_PANELS
        || count.checked_mul(4_usize.pow(levels)).is_none_or(|n|n>maximum)
        || boundary.weights.len()!=boundary.state_modes.len() || boundary.weights.is_empty()
        || boundary.weights.iter().any(|r|r.len()!=count || r.iter().any(|v|!v.is_finite()))
        || boundary.triangles.iter().flatten().flatten().any(|v|!v.is_finite()) {
        return Err("radiation refinement exceeds its complete geometry/source budget".into());
    }
    let mut points=Vec::new(); let mut lookup=BTreeMap::new(); let mut faces=Vec::with_capacity(count);
    for (parent,t) in boundary.triangles.iter().enumerate() {
        if !area(*t).is_finite() || area(*t)<=0.0 { return Err("degenerate radiation panel".into()); }
        let nodes=t.map(|p|*lookup.entry(p.map(|x|if x==0.0 {0}else{x.to_bits()})).or_insert_with(|| {
            let i=points.len();points.push(p);i
        }));
        faces.push(Face {nodes,chart:None,parent});
    }
    check_closed(&faces.iter().map(|f|f.nodes).collect::<Vec<_>>())?;
    for _ in 0..levels {
        if gate.is_requested() { return Err("radiation refinement cancelled".into()); }
        let mut marked=BTreeSet::new();
        for f in &faces { for i in 0..3 { marked.insert(edge(f.nodes[i],f.nodes[(i+1)%3])); } }
        faces=refine(&mut points,&faces,&marked,maximum)?;
    }
    check_closed(&faces.iter().map(|f|f.nodes).collect::<Vec<_>>())?;
    Ok(Boundary {triangles:faces.iter().map(|f|f.nodes.map(|i|points[i])).collect(),
        weights:boundary.weights.iter().map(|r|faces.iter().map(|f|r[f.parent]).collect()).collect(),
        state_modes:boundary.state_modes.clone()})
}
