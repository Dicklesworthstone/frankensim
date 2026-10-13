//! P1 pressure geometry and analytic Laplace moments, not a second mesh owner.
use super::*;
use std::collections::BTreeMap;

pub(super) type Point = [f64; 3];
pub(super) fn sub(a: Point, b: Point) -> Point { core::array::from_fn(|i| a[i]-b[i]) }
pub(super) fn dot(a: Point, b: Point) -> f64 { a.iter().zip(b).map(|(x,y)| x*y).sum() }
pub(super) fn cross(a: Point, b: Point) -> Point {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
pub(super) fn norm(a: Point) -> f64 { a[0].hypot(a[1]).hypot(a[2]) }
pub(super) fn add_scaled(a: &mut Point, b: Point, scale: f64) {
    for i in 0..3 { a[i] += b[i]*scale; }
}
pub(super) const X4: [f64;4] = [0.06943184420297371,0.33000947820757187,0.6699905217924281,0.9305681557970262];
pub(super) const W4: [f64;4] = [0.17392742256872692,0.32607257743127307,0.32607257743127307,0.17392742256872692];
pub(super) const X8: [f64;8] = [0.019855071751231884,0.10166676129318664,0.2372337950418355,0.4082826787521751,
    0.5917173212478249,0.7627662049581645,0.8983332387068134,0.9801449282487681];
pub(super) const W8: [f64;8] = [0.05061426814518813,0.11119051722668724,0.15685332293894366,0.181341891689181,
    0.181341891689181,0.15685332293894366,0.11119051722668724,0.05061426814518813];

#[derive(Clone, Copy)]
pub(super) struct Sample { pub point: Point, pub shape: [f64;3], pub weight: f64 }
pub(super) struct Facet {
    pub triangle: [Point;3], pub ids: [usize;3], pub normal: Point, pub area: f64,
    pub gradient: [Point;3], pub curl: [Point;3], pub samples: [Sample;16],
}
impl Facet {
    fn new(triangle: [Point;3], ids: [usize;3], area:f64, normal:Point) -> Result<Self,HelmholtzError> {
        // Preserve the geometry owner's exact area/normal bytes: the caller's
        // work projection uses these same areas, not a separately rounded norm.
        let size=2.0*area;
        if !size.is_finite() || area<=0.0 || normal.iter().any(|v|!v.is_finite()) {
            return Err(bad("degenerate reciprocal boundary facet"));
        }
        let gradient=core::array::from_fn(|i|cross(normal,
            sub(triangle[(i+2)%3],triangle[(i+1)%3])).map(|v|v/size));
        let curl=gradient.map(|g|cross(normal,g));
        let samples=core::array::from_fn(|i| {
            let u=X4[i/4]; let v=X4[i%4]; let shape=[(1.0-u)*(1.0-v),u,(1.0-u)*v];
            let point=core::array::from_fn(|d| (0..3).map(|a|shape[a]*triangle[a][d]).sum());
            Sample {point,shape,weight:size*(1.0-u)*W4[i/4]*W4[i%4]}
        });
        if gradient.iter().flatten().any(|v|!v.is_finite()) { return Err(bad("unrepresentable reciprocal boundary gradient")); }
        Ok(Self {triangle,ids,normal,area,gradient,curl,samples})
    }
}

pub(super) fn facets(surface:&SpherePanels) -> Result<(Vec<Facet>,usize,Point,f64),HelmholtzError> {
    let triangles=surface.triangles().ok_or_else(||bad("reciprocal radiation requires retained triangles"))?;
    let mut nodes=BTreeMap::<[u64;3],usize>::new();
    let mut edges=BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
    let mut lower=[f64::INFINITY;3]; let mut upper=[f64::NEG_INFINITY;3];
    let mut out=Vec::with_capacity(triangles.len());
    for (index,&triangle) in triangles.iter().enumerate() {
        if triangle.iter().flatten().any(|v|!v.is_finite() || v.abs()>1e9) {return Err(bad("nonfinite reciprocal boundary"));}
        let ids=triangle.map(|p| {
            for d in 0..3 {lower[d]=lower[d].min(p[d]);upper[d]=upper[d].max(p[d]);}
            let key=p.map(|v|if v==0.0 {0}else{v.to_bits()}); let next=nodes.len();
            *nodes.entry(key).or_insert(next)
        });
        for a in 0..3 {let i=ids[a];let j=ids[(a+1)%3];edges.entry((i.min(j),i.max(j))).or_default().push((i,j));}
        out.push(Facet::new(triangle,ids,surface.areas()[index],surface.normals()[index])?);
    }
    if out.len()<4 || edges.values().any(|e|e.len()!=2 || e[0]!=(e[1].1,e[1].0)) {
        return Err(bad("reciprocal radiation requires an oriented closed conforming surface"));
    }
    let widths=core::array::from_fn(|d| (upper[d]-lower[d]).next_up());
    let diameter=norm(widths);
    if !diameter.is_finite() || diameter<=0.0 {return Err(bad("reciprocal boundary extent"));}
    Ok((out,nodes.len(),widths,diameter))
}

/// Static single-layer integral and its three P1 moments, and P1 double-layer
/// moments. Extend the existing constant-triangle edge/solid-angle formula:
/// integral (y-x_parallel)/r = contour integral r nu, and
/// integral (y-x_parallel)*h/r^3 = -h * contour integral nu/r.
/// These keep self and nearly touching interactions out of area sampling.
pub(super) fn static_moments(x:Point, f:&Facet, same:bool) -> (f64,[f64;3],[f64;3]) {
    let relative=f.triangle.map(|p|sub(p,x)); let radii=relative.map(norm);
    let h=if same {0.0}else{-dot(relative[0],f.normal)};
    let den=radii[0]*radii[1]*radii[2]+dot(relative[0],relative[1])*radii[2]
        +dot(relative[1],relative[2])*radii[0]+dot(relative[2],relative[0])*radii[1];
    let omega=if h==0.0 {0.0}else{2.0*det::atan2(h*2.0*f.area,den)};
    let mut edge_sum=0.0; let mut moment=[0.0;3]; let mut derivative=[0.0;3];
    for a in 0..3 {
        let b=(a+1)%3;let d=sub(f.triangle[b],f.triangle[a]);let length=norm(d);
        let tangent=d.map(|v|v/length);let outward=cross(tangent,f.normal);
        let start=dot(relative[a],tangent);let end=start+length;
        let perpendicular=norm(cross(tangent,relative[a]));
        let line=if start>=0.0 {det::ln((radii[b]+end)/(radii[a]+start))}
            else if end<=0.0 {det::ln((radii[a]-start)/(radii[b]-end))}
            else {det::ln((radii[a]-start)/perpendicular)+det::ln((radii[b]+end)/perpendicular)};
        edge_sum+=dot(relative[a],outward)*line;
        add_scaled(&mut derivative,outward,line);
        add_scaled(&mut moment,outward,0.5*(end*radii[b]-start*radii[a]+perpendicular*perpendicular*line));
    }
    let single=(edge_sum-h*omega)/(4.0*PI);
    let mut phi=f.gradient.map(|g|dot(g,sub(x,f.triangle[0])));phi[0]+=1.0;
    let s=core::array::from_fn(|a|phi[a]*single+dot(f.gradient[a],moment)/(4.0*PI));
    let d=core::array::from_fn(|a|phi[a]*omega/(4.0*PI)-h*dot(f.gradient[a],derivative)/(4.0*PI));
    (single,s,d)
}

#[derive(Clone, Copy)]
pub(super) struct Pair { pub single:f64, pub shape:[[f64;3];3] }
fn directed(a:&Facet,b:&Facet,same:bool) -> (Pair,[f64;3]) {
    let mut result=Pair {single:0.0,shape:[[0.0;3];3]};let mut k=[0.0;3];
    for x in a.samples {
        let (s,u,d)=static_moments(x.point,b,same);
        result.single+=x.weight*s;
        for i in 0..3 {
            k[i]+=x.weight*d[i];
            for j in 0..3 {result.shape[i][j]+=x.weight*x.shape[i]*u[j];}
        }
    }
    (result,k)
}
pub(super) fn pair(a:&Facet,b:&Facet,same:bool) -> (Pair,[f64;3],[f64;3]) {
    let (forward,ka)=directed(a,b,same);
    let (back,kb)=if same {(forward,ka)}else{directed(b,a,false)};
    // Symmetric testing/source quadrature of the reciprocal integral, not a
    // repair to a computed impedance or an alteration of a negative eigenvalue.
    (Pair {single:0.5*forward.single+0.5*back.single,
        shape:core::array::from_fn(|i|core::array::from_fn(|j|0.5*forward.shape[i][j]+0.5*back.shape[j][i]))},ka,kb)
}
