//! A prescribed volume-flow aperture on the existing closed exterior BEM mesh.
//! This is a Neumann source boundary, not radiation feedback on the neck.
use super::{Boundary, Error, MAX_INPUTS, MAX_PANELS, check_closed};
use fs_couple::render::plate::impact::cavity::neck::NeckRadiationPort;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

type Point = [f64; 3];
type Key = [u64; 3];
fn key(p: Point) -> Key { p.map(|x| if x == 0.0 { 0 } else { x.to_bits() }) }
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 { a[0]*b[1]-a[1]*b[0] }
fn sub(a: Point,b: Point)->Point {std::array::from_fn(|i|a[i]-b[i])}
fn dot(a: Point,b: Point)->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn area(t: [Point; 3])->f64 {
    let a=sub(t[1],t[0]);let b=sub(t[2],t[0]);
    let n=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
    0.5*dot(n,n).sqrt()
}

/// Explicit model selection. A bare vent/audio request still refuses.
pub fn option(args:&mut Vec<String>)->Result<bool,Error> {
    let n=args.iter().filter(|a|a.as_str()=="--prescribed-vent-radiation").count();
    if n>1 {return Err("--prescribed-vent-radiation may be supplied only once".into());}
    args.retain(|a|a!="--prescribed-vent-radiation");Ok(n==1)
}

// Identify a planar vertical wall strip by its two exact XY endpoint addresses.
fn wall(t: &[Point; 3])->Option<([f64; 2],[f64; 2])> {
    let a=[t[0][0],t[0][1]];
    let b=t.iter().map(|p|[p[0],p[1]]).find(|p|*p!=a)?;
    if t.iter().any(|p|[p[0],p[1]]!=a && [p[0],p[1]]!=b) {return None;}
    Some(if a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])).is_lt() {(a,b)}else{(b,a)})
}

// Preserve each outer boundary vertex byte-for-byte. Adjacent wall/head panels
// must not acquire T-junctions or cracks when the local strip is retriangulated.
fn perimeter(triangles:&[[Point; 3]])->Result<Vec<Point>,Error> {
    let mut edges=BTreeMap::<(Key,Key),Vec<(Point,Point)>>::new();
    for t in triangles {for i in 0..3 {
        let (a,b)=(t[i],t[(i+1)%3]);let (ka,kb)=(key(a),key(b));
        edges.entry((ka.min(kb),ka.max(kb))).or_default().push((a,b));
    }}
    let mut next=BTreeMap::<Key,(Point,Point)>::new();
    for e in edges.values() {match e.as_slice() {
        &[(a,b)]=>{if next.insert(key(a),(a,b)).is_some() {return Err("branching aperture perimeter".into());}},
        &[(a,b),(c,d)] if key(a)==key(d) && key(b)==key(c)=>{},
        _=>return Err("aperture strip has inconsistent surface edges".into()),
    }}
    let start=*next.keys().next().ok_or("empty aperture perimeter")?;
    let mut at=start;let mut out=Vec::with_capacity(next.len());
    for _ in 0..next.len() {
        let &(a,b)=next.get(&at).ok_or("open aperture perimeter")?;out.push(a);at=key(b);
        if at==start && out.len()!=next.len() {return Err("multiple aperture perimeter loops".into());}
    }
    if at!=start {return Err("unclosed aperture perimeter".into());}Ok(out)
}

struct Patch { triangles:Vec<[Point; 3]>, first_open:usize, area_m2:f64 }

// A convex planar wall polygon minus an inscribed circular polygon, followed by
// the aperture disk. Radial fans retain every collinear outer edge subdivision.
// The disk's volume integral, not peak velocity, fixes its discrete source scale.
fn partition(outer:&[Point],center:Point,tangent:Point,vertical:Point,radius:f64,segments:usize)
    ->Result<Patch,Error>
{
    if outer.len()<3 || !(12..=128).contains(&segments) {return Err("aperture partition budget".into());}
    let local:Vec<_>=outer.iter().map(|&p|{let d=sub(p,center);[dot(d,tangent),dot(d,vertical)]}).collect();
    for i in 0..local.len() {
        let (a,b)=(local[i],local[(i+1)%local.len()]);let edge=[b[0]-a[0],b[1]-a[1]];
        let distance=cross(edge,[-a[0],-a[1]])/edge[0].hypot(edge[1]);
        if !distance.is_finite() || distance<=radius*(1.0+1e-10) {
            return Err("circular aperture crosses a faceted wall seam or head/rim; its supplied position is never snapped".into());
        }
    }
    let outer_angles:Vec<_>=local.iter().map(|p|p[1].atan2(p[0]).rem_euclid(TAU)).collect();
    // Each annular sector is concave on its circular side. A single fan from
    // one outer corner can cross the aperture. Split at the edge-normal ray,
    // which both endpoints see because the edge is strictly outside the disk.
    let split_angles:Vec<_>=(0..local.len()).map(|i| {
        let (a,b)=(local[i],local[(i+1)%local.len()]);
        let normal=[b[1]-a[1],a[0]-b[0]];
        let relative=cross(a,normal).atan2(a[0]*normal[0]+a[1]*normal[1]);
        let span=(outer_angles[(i+1)%local.len()]-outer_angles[i]).rem_euclid(TAU);
        (outer_angles[i]+relative.clamp(0.0,span)).rem_euclid(TAU)
    }).collect();
    let mut angles:Vec<_>=(0..segments).map(|i|TAU*i as f64/segments as f64)
        .chain(outer_angles.iter().copied()).chain(split_angles.iter().copied()).collect();
    angles.sort_by(f64::total_cmp);angles.dedup_by(|a,b|(*a-*b).abs()<1e-12);
    if angles.len()>1 && TAU-angles[angles.len()-1]+angles[0]<1e-12 {angles.pop();}
    let ring:Vec<Point>=angles.iter().map(|a|std::array::from_fn(|i|
        center[i]+radius*(a.cos()*tangent[i]+a.sin()*vertical[i]))).collect();
    let address=|a:&f64|angles.iter().position(|b| {
        let d=(a-b).abs();d.min(TAU-d)<1e-12
    }).ok_or("missing aperture radial address");
    let indices=outer_angles.iter().map(address).collect::<Result<Vec<_>,_>>()?;
    let split_indices=split_angles.iter().map(address).collect::<Result<Vec<_>,_>>()?;
    let mut triangles=Vec::new();
    for i in 0..outer.len() {
        let next=(i+1)%outer.len();let split=split_indices[i];
        triangles.push([outer[i],outer[next],ring[split]]);
        for (anchor,start,end) in [(outer[i],indices[i],split),(outer[next],split,indices[next])] {
            let mut at=start;let mut walked=0;
            while at!=end {
                let after=(at+1)%ring.len();triangles.push([anchor,ring[after],ring[at]]);
                at=after;walked+=1;if walked>ring.len() {return Err("invalid aperture angular order".into());}
            }
        }
    }
    let first_open=triangles.len();
    for i in 0..ring.len() {triangles.push([center,ring[i],ring[(i+1)%ring.len()]]);}
    if triangles.iter().any(|t|!area(*t).is_finite() || area(*t)<=0.0) {
        return Err("degenerate aperture panel".into());
    }
    let area_m2=triangles[first_open..].iter().map(|t|area(*t)).sum();
    Ok(Patch {triangles,first_open,area_m2})
}

impl Boundary {
    /// Append one physical neck source, preserving the head-mode prefix and all
    /// prior normal-velocity weights. The patch closes the computational surface
    /// with prescribed OUTWARD flow; it is not a rigid cap and not a monopole mix.
    /// The circular aperture must lie within one actual planar sidewall facet.
    /// It may cross that facet's axial subdivisions, which are retriangulated.
    /// Failure leaves no partially usable candidate. No mechanical state changes.
    pub(crate) fn with_sidewall_aperture(self,port:NeckRadiationPort,azimuth:f64,
        axial_from_top:f64,depth:f64)->Result<Self,Error>
    {
        let radius=(port.area_m2/PI).sqrt();let medium=fs_bem::helmholtz::Medium::air();
        if ![radius,port.volume_weight_m2_per_sqrt_kg,port.effective_length_m,depth].iter().all(|x|x.is_finite() && *x>0.0)
            || !azimuth.is_finite() || !axial_from_top.is_finite()
            || self.weights.len()>=MAX_INPUTS || self.state_modes.contains(&port.coordinate)
            || self.weights.len()!=self.state_modes.len() || self.triangles.len()>MAX_PANELS
            || self.weights.iter().any(|r|r.len()!=self.triangles.len()) {
            return Err("prescribed aperture needs a finite admitted neck and unused source capacity/address".into());
        }
        // The uniform-velocity slug must cover the SAME acoustic bake window,
        // not just the lower cavity-mode window used by the mechanical builder.
        if TAU*1640.0/medium.sound_speed*radius.max(port.effective_length_m)>0.3 {
            return Err("neck exceeds the compact chart over the 40..1640 Hz radiation window".into());
        }
        let direction=[azimuth.cos(),azimuth.sin()];let mut selected=None;
        for t in &self.triangles {if let Some((a,b))=wall(t) {
            let d=[b[0]-a[0],b[1]-a[1]];let denominator=cross(d,direction);
            if denominator==0.0 {continue;}
            let fraction=-cross(a,direction)/denominator;
            let p=[a[0]+fraction*d[0],a[1]+fraction*d[1]];
            if fraction>0.0 && fraction<1.0 && p[0]*direction[0]+p[1]*direction[1]>0.0 {
                selected=Some((a,b,p));break;
            }
        }}
        let (a,b,p)=selected.ok_or("aperture azimuth has no interior planar sidewall facet")?;
        let picked:Vec<_>=self.triangles.iter().enumerate().filter_map(|(i,t)|
            (wall(t)==Some((a,b))).then_some(i)).collect();
        if picked.is_empty() || self.weights.iter().any(|r|picked.iter().any(|&i|r[i]!=0.0)) {
            return Err("aperture must replace a rigid sidewall, not a moving head or shell mode".into());
        }
        let outer=perimeter(&picked.iter().map(|&i|self.triangles[i]).collect::<Vec<_>>())?;
        let length=(b[0]-a[0]).hypot(b[1]-a[1]);let tangent=[(b[0]-a[0])/length,(b[1]-a[1])/length,0.0];
        // tangent cross vertical must point outward.
        let sign=if tangent[1]*p[0]-tangent[0]*p[1]>0.0 {1.0}else{-1.0};
        let patch=partition(&outer,[p[0],p[1],0.5*depth-axial_from_top],tangent,[0.0,0.0,sign],radius,32)?;
        if (patch.area_m2-port.area_m2).abs()>0.01*port.area_m2 {
            return Err("aperture polygon exceeds its one-percent geometric area budget".into());
        }
        let count=self.triangles.len()-picked.len()+patch.triangles.len();
        if count>MAX_PANELS {return Err("resolved aperture exceeds the existing BEM panel budget".into());}
        let mut triangles=Vec::with_capacity(count);let mut weights=vec![Vec::with_capacity(count);self.weights.len()];
        for (i,t) in self.triangles.into_iter().enumerate() {if picked.binary_search(&i).is_err() {
            triangles.push(t);for (out,row) in weights.iter_mut().zip(&self.weights) {out.push(row[i]);}
        }}
        let first_open=triangles.len()+patch.first_open;triangles.extend(patch.triangles);
        for row in &mut weights {row.resize(count,0.0);}
        let gain=port.volume_weight_m2_per_sqrt_kg/patch.area_m2;
        if !gain.is_finite() || gain<=0.0 {return Err("unrepresentable aperture volume-flow projection".into());}
        let mut flow=vec![0.0;count];flow[first_open..].fill(gain);weights.push(flow);
        let mut ids=BTreeMap::<Key,usize>::new();
        let indices:Vec<_>=triangles.iter().map(|t|t.map(|p| {let next=ids.len();*ids.entry(key(p)).or_insert(next)})).collect();
        check_closed(&indices)?;
        let mut state_modes=self.state_modes;state_modes.push(port.coordinate);
        eprintln!("prescribed-flow aperture: nominal_area_m2={}, polygon_area_m2={}, added_source_coordinate={}; existing head sources and actual neck volume flow, no radiation backreaction or added end correction",port.area_m2,patch.area_m2,port.coordinate);
        Ok(Self {triangles,weights,state_modes})
    }
}

#[cfg(test)]
#[path = "aperture/tests.rs"]
mod tests;
