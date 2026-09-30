//! A usable source-derived Model D soundboard, not an empty measurement slot.
//!
//! Geometry: Boutillon, Ege & Paulello, Acoustics 2012, Fig. 2 (Claire Pichet's
//! technical drawing), https://arxiv.org/abs/1210.3948 . Coordinates below are
//! our approximate manual transcription of the geometric facts in that figure;
//! the copyrighted drawing is NOT redistributed. Table 1's 0.127 m mean rib
//! spacing sets the scale. Ribs are made parallel in this reduced reconstruction.
//! The 17 heights follow the diagram's H labels; rib widths and end tapers,
//! bridge widths, cut-off bar section, damping, G and nu are explicit estimates.
//!
//! Materials/dimensions: https://www.steinway.com/pianos/steinway/grand/model-d
//! Sitka spruce panel (9 -> 6 mm), sugar pine ribs, hard-maple bridges. Panel
//! elastic constants use the paper's *average spruce*, not a specimen fit.
//! Grain follows the long bridge, perpendicular to the ribs. Crown, rim/plate
//! flexibility, glue slip and downbearing remain outside this flat-plate image.
//!
//! Produces BOTH the existing fs-plate input and an OBJ from the SAME vertices,
//! sections and beams. All work is cold. No eigensolver or contact law lives here.
use std::collections::BTreeMap;
use std::fmt::Write;

const SLOPE: f64 = 0.8;
const FIRST_RIB: f64 = 241.0;
const LAST_RIB: f64 = 1439.0;
const SPACING_M: f64 = 0.127;
// Coordinates in the 722 x 1000 raster representation of Fig. 2, not metres.
const OUTLINE: &[[f64; 2]] = &[
    [57.,100.],[61.,75.],[89.,54.],[137.,44.],[186.,41.],
    [243.,42.],[285.,50.],[334.,65.],[378.,90.],[399.,129.],
    [410.,220.],[420.,306.],[432.,391.],[445.,461.],[456.,540.],
    [467.,616.],[487.,671.],[508.,710.],[550.,757.],[590.,789.],
    [638.,825.],[676.,861.],[677.,970.],[55.,920.],[57.,730.],
];
// (transverse sweep station, maximum height in mm), ribs 17 down to 1.
const RIBS: &[(f64, f64)] = &[
    (241.,27.),(328.,27.),(402.,26.),(481.,29.),(565.,26.),
    (648.,26.),(732.,23.),(816.,26.),(891.,25.),(969.,26.),
    (1048.,23.),(1118.,24.),(1187.,23.),(1254.,22.),
    (1321.,19.),(1387.,19.),(1439.,17.),
];
// Bass end -> treble end. These describe a bridge CURVE, not the string scale.
const LONG_BRIDGE: &[[f64; 2]] = &[
    [150.,153.],[132.,178.],[128.,212.],[138.,252.],
    [156.,300.],[178.,350.],[197.,401.],[217.,450.],[242.,500.],
    [259.,550.],[287.,607.],[320.,663.],[352.,718.],[393.,772.],
    [431.,825.],[472.,868.],[519.,905.],[570.,932.],[624.,953.],[670.,960.],
];
const BASS_BRIDGE: &[[f64; 2]] = &[
    [150.,153.],[184.,147.],[224.,165.],[260.,193.],[286.,239.],[302.,294.],
    [315.,358.],[326.,429.],
];
const CUTOFF: &[[f64; 2]] = &[[57.,730.],[193.,930.]];

// RT-0425 Appendix A, C1..B7 (MIDI 24..107): bridge coupling centre (x0,y0)
// in metres. These belong to its reference board, not to this reconstructed
// outline; project onto the nearest point of the corresponding structural
// bridge below. A0..B0 and C8 are extrapolated because Appendix A omits them.
const RT0425_BRIDGE: [[f64; 2]; 84] = [
    [0.40,1.69], [0.42,1.67], [0.44,1.65], [0.46,1.63], [0.48,1.60], [0.50,1.56],
    [0.52,1.52], [0.53,1.48], [0.55,1.44], [0.56,1.39], [0.57,1.35], [0.58,1.30],
    [0.59,1.26], [0.60,1.22], [0.61,1.18], [0.61,1.14], [0.62,1.11], [0.23,1.58],
    [0.25,1.49], [0.27,1.41], [0.30,1.34], [0.32,1.26], [0.35,1.19], [0.37,1.13],
    [0.40,1.06], [0.42,1.01], [0.45,0.95], [0.47,0.90], [0.49,0.85], [0.51,0.80],
    [0.54,0.75], [0.56,0.71], [0.58,0.67], [0.60,0.63], [0.62,0.60], [0.64,0.56],
    [0.66,0.53], [0.68,0.50], [0.70,0.47], [0.71,0.44], [0.73,0.41], [0.75,0.39],
    [0.76,0.36], [0.78,0.34], [0.80,0.32], [0.81,0.30], [0.83,0.28], [0.85,0.27],
    [0.86,0.25], [0.88,0.23], [0.90,0.22], [0.91,0.21], [0.93,0.19], [0.94,0.18],
    [0.96,0.17], [0.97,0.16], [0.99,0.15], [1.01,0.14], [1.02,0.13], [1.04,0.12],
    [1.05,0.11], [1.07,0.10], [1.09,0.09], [1.10,0.08], [1.12,0.07], [1.13,0.06],
    [1.15,0.06], [1.16,0.05], [1.18,0.04], [1.19,0.04], [1.21,0.03], [1.22,0.02],
    [1.24,0.02], [1.25,0.01], [1.27,0.01], [1.28,0.00], [1.29,0.00], [1.31,-0.01],
    [1.32,-0.01], [1.33,-0.02], [1.34,-0.02], [1.36,-0.02], [1.37,-0.03], [1.38,-0.03],
];

pub struct Preset {
    pub geometry: String,
    pub obj: String,
}
#[derive(Clone)]
struct Row { station: f64, nodes: Vec<usize>, features: [Option<usize>; 3] }
struct Beam { a: usize, b: usize, width: f64, height: f64, z: f64,
    e: f64, g: f64, rho: f64, group: String }

fn station(p: [f64; 2]) -> f64 { p[1] + SLOPE * p[0] }
fn scale() -> f64 { SPACING_M * (1.0 + SLOPE*SLOPE).sqrt() * 16.0 / (LAST_RIB-FIRST_RIB) }
fn si(p: [f64; 2]) -> [f64; 2] { [(p[0]-55.0)*scale(), (970.0-p[1])*scale()] }
fn area2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn intersections(curve: &[[f64; 2]], t: f64, closed: bool) -> Vec<f64> {
    let mut xs = Vec::new();
    let count = if closed { curve.len() } else { curve.len()-1 };
    for i in 0..count {
        let a=curve[i]; let b=curve[(i+1)%curve.len()];
        let (ta,tb)=(station(a),station(b));
        if (tb-ta).abs()<1e-10 { continue; }
        let u=(t-ta)/(tb-ta);
        if (-1e-10..=1.0+1e-10).contains(&u) { xs.push(a[0]+u.clamp(0.0,1.0)*(b[0]-a[0])); }
    }
    xs.sort_by(f64::total_cmp); xs.dedup_by(|a,b| (*a-*b).abs()<1e-8); xs
}
fn distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let d=[b[0]-a[0],b[1]-a[1]];
    let u=(((p[0]-a[0])*d[0]+(p[1]-a[1])*d[1])/(d[0]*d[0]+d[1]*d[1])).clamp(0.,1.);
    ((p[0]-a[0]-u*d[0]).powi(2)+(p[1]-a[1]-u*d[1]).powi(2)).sqrt()
}
fn thickness(p: [f64; 2]) -> f64 {
    let d=(0..OUTLINE.len()).map(|i| distance(p,si(OUTLINE[i]),si(OUTLINE[(i+1)%OUTLINE.len()])))
        .fold(f64::INFINITY,f64::min);
    // Authored smooth taper profile; only its 9/6 mm endpoints are published.
    0.006+0.003*(d/0.20).clamp(0.0,1.0)
}
fn append_triangle(tris: &mut Vec<[usize; 3]>, xy: &[[f64; 2]], a:usize,b:usize,c:usize) {
    if area2(xy[a],xy[b],xy[c])>0.0 { tris.push([a,b,c]); } else { tris.push([a,c,b]); }
}
fn curve_point(curve: &[[f64;2]], fraction:f64)->[f64;2] {
    let lengths:Vec<f64>=curve.windows(2).map(|p| ((p[1][0]-p[0][0]).powi(2)+(p[1][1]-p[0][1]).powi(2)).sqrt()).collect();
    let mut remaining=fraction*lengths.iter().sum::<f64>();
    for (p,len) in curve.windows(2).zip(lengths) {
        if remaining<=len { return [p[0][0]+remaining/len*(p[1][0]-p[0][0]),p[0][1]+remaining/len*(p[1][1]-p[0][1])]; }
        remaining-=len;
    }
    *curve.last().expect("nonempty source curve")
}
fn nearest_bridge_point(p: [f64; 2], curve: &[[f64; 2]]) -> ([f64; 2], f64) {
    let mut best = ([0.0, 0.0], f64::INFINITY);
    for segment in curve.windows(2) {
        let (a, b) = (si(segment[0]), si(segment[1]));
        let d = [b[0] - a[0], b[1] - a[1]];
        let t = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1])
            / (d[0] * d[0] + d[1] * d[1])).clamp(0.0, 1.0);
        let q = [a[0] + t * d[0], a[1] + t * d[1]];
        let error = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
        if error < best.1 { best = (q, error); }
    }
    best
}
fn published_bridge_point(key: u8) -> ([f64; 2], f64) {
    let p = match key {
        21..=23 => {
            let steps = f64::from(24 - key);
            [RT0425_BRIDGE[0][0] - steps * 0.02,
             RT0425_BRIDGE[0][1] + steps * 0.02]
        }
        24..=107 => RT0425_BRIDGE[usize::from(key - 24)],
        108 => [1.39, -0.04],
        _ => unreachable!("piano key range is A0..C8"),
    };
    nearest_bridge_point(p, if key <= 40 { BASS_BRIDGE } else { LONG_BRIDGE })
}
fn locate(p:[f64;2],xy:&[[f64;2]],tris:&[[usize;3]])->Result<(usize,[f64;3]),String> {
    for (i,&[a,b,c]) in tris.iter().enumerate() {
        let twice=area2(xy[a],xy[b],xy[c]);
        let mut w=[area2(p,xy[b],xy[c])/twice,area2(xy[a],p,xy[c])/twice,area2(xy[a],xy[b],p)/twice];
        if w.iter().all(|x| *x>=-1e-9 && *x<=1.0+1e-9) {
            for x in &mut w { *x=x.clamp(0.0,1.0); }
            let sum=w.iter().sum::<f64>(); for x in &mut w { *x/=sum; }
            return Ok((i,w));
        }
    }
    Err("source bridge station falls outside reconstructed panel".into())
}

/// Rib-aligned, conforming sweep of this particular drawing. Feature points are
/// inserted in each row, so structural bridges/cut-off and ribs are actual mesh
/// beams, not nearest-node decorations. `divisions` controls transverse refinement.
pub fn build(divisions: usize) -> Result<Preset,String> {
    build_with_contacts(divisions, false)
}

pub fn build_with_rt0425_contacts(divisions: usize) -> Result<Preset,String> {
    build_with_contacts(divisions, true)
}

fn build_with_contacts(divisions: usize, published_contacts: bool) -> Result<Preset,String> {
    if !(4..=32).contains(&divisions) { return Err("Model D mesh divisions must be 4..32".into()); }
    build_inner(divisions, published_contacts, None)
}

#[cfg(test)]
pub(crate) fn build_probe(divisions: usize) -> Result<Preset,String> {
    if !(33..=64).contains(&divisions) { return Err("research mesh divisions must be 33..64".into()); }
    build_inner(divisions, false, None)
}

#[cfg(test)]
pub(crate) fn build_probe_refined(divisions: usize, maximum_station_gap: f64) -> Result<Preset,String> {
    if !(33..=64).contains(&divisions) || !maximum_station_gap.is_finite()
        || maximum_station_gap<=0. {
        return Err("invalid research mesh divisions or row gap".into());
    }
    build_inner(divisions, false, Some(maximum_station_gap))
}

fn build_inner(divisions: usize, published_contacts: bool, maximum_station_gap: Option<f64>) -> Result<Preset,String> {
    let curves=[LONG_BRIDGE,BASS_BRIDGE,CUTOFF];
    let mut levels:Vec<f64>=OUTLINE.iter().copied().map(station).collect();
    levels.extend(RIBS.iter().map(|r|r.0));
    for curve in curves { levels.extend(curve.iter().copied().map(station)); }
    levels.sort_by(f64::total_cmp); levels.dedup_by(|a,b|(*a-*b).abs()<1e-8);
    let mut extra=Vec::new();
    for p in levels.windows(2) { if p[1]-p[0]>55.0 {extra.push((p[0]+p[1])*0.5);} }
    levels.extend(extra); levels.sort_by(f64::total_cmp);
    if let Some(maximum_gap)=maximum_station_gap {
        let mut refined=Vec::new();
        for pair in levels.windows(2) {
            let segments=((pair[1]-pair[0])/maximum_gap).ceil() as usize;
            for index in 1..segments {
                refined.push(pair[0]+(pair[1]-pair[0])*index as f64/segments as f64);
            }
        }
        levels.extend(refined); levels.sort_by(f64::total_cmp);
    }
    let mut xy=Vec::new(); let mut rows=Vec::new();
    for t in levels {
        let ends=intersections(OUTLINE,t,true);
        if ends.is_empty()||ends.len()>2 {return Err("source outline is not sweep-monotone".into());}
        let lo=ends[0]; let hi=*ends.last().expect("nonempty span");
        let mut xs=if hi-lo<1e-8 {vec![lo]} else {(0..=divisions).map(|i|lo+(hi-lo)*i as f64/divisions as f64).collect()};
        let mut feature_x=[None;3];
        for (k,curve) in curves.iter().enumerate() {
            let hit=intersections(curve,t,false);
            if hit.len()>1 {return Err("source feature is not sweep-monotone".into());}
            if let Some(&x)=hit.first() {
                if x<lo-1e-7||x>hi+1e-7 {return Err("source feature lies outside panel".into());}
                xs.push(x.clamp(lo,hi));feature_x[k]=Some(x.clamp(lo,hi));
            }
        }
        xs.sort_by(f64::total_cmp);xs.dedup_by(|a,b|(*a-*b).abs()<1e-8);
        let mut row=Row{station:t,nodes:Vec::new(),features:[None;3]};
        for x in xs {
            let id=xy.len();xy.push(si([x,t-SLOPE*x]));row.nodes.push(id);
            for k in 0..3 {if feature_x[k].is_some_and(|p|(p-x).abs()<1e-8){row.features[k]=Some(id);}}
        }
        rows.push(row);
    }
    let mut tris=Vec::new();
    // Each bridge/cut-off beam between rows must also be a triangle edge.
    // Partition the strip at features shared by both rows, then merge the
    // normalized abscissae within each partition. Merging entire rows lets
    // the chosen diagonal cross a structural beam, disconnecting its span
    // from the plate interior except at the two endpoint nodes.
    for pair in rows.windows(2) {
        let (a,b)=(&pair[0].nodes,&pair[1].nodes);
        let mut cuts=vec![(0,0)];
        for k in 0..3 {
            if let (Some(left),Some(right))=(pair[0].features[k],pair[1].features[k]) {
                let ia=a.iter().position(|&node|node==left).ok_or("missing upper feature node")?;
                let ib=b.iter().position(|&node|node==right).ok_or("missing lower feature node")?;
                cuts.push((ia,ib));
            }
        }
        cuts.push((a.len()-1,b.len()-1));
        cuts.sort_unstable();cuts.dedup();
        for limits in cuts.windows(2) {
            let ((a0,b0),(a1,b1))=(limits[0],limits[1]);
            if b0>b1 {return Err("structural feature edges cross inside a sweep strip".into());}
            let (left,right)=(&a[a0..=a1],&b[b0..=b1]);
            let (mut i,mut j)=(0,0);
            while i+1<left.len()||j+1<right.len() {
                let next=|ids:&[usize],k:usize| {
                    if k+1==ids.len(){f64::INFINITY}else{
                        (xy[ids[k+1]][0]-xy[ids[0]][0])/(xy[*ids.last().unwrap()][0]-xy[ids[0]][0])
                    }
                };
                if next(left,i)<=next(right,j) {
                    append_triangle(&mut tris,&xy,left[i],left[i+1],right[j]);i+=1;
                } else {
                    append_triangle(&mut tris,&xy,left[i],right[j+1],right[j]);j+=1;
                }
            }
        }
    }
    let mut edges:BTreeMap<(usize,usize),usize>=BTreeMap::new();
    for &[a,b,c] in &tris {
        if area2(xy[a],xy[b],xy[c])<=1e-15 {return Err("degenerate reconstruction triangle".into());}
        for (u,v) in [(a,b),(b,c),(c,a)] {*edges.entry((u.min(v),u.max(v))).or_default()+=1;}
    }
    let mut beams=Vec::new();
    for (index,&(t,max_h)) in RIBS.iter().enumerate() {
        let row=rows.iter().find(|r|(r.station-t).abs()<1e-8).ok_or("missing rib row")?;
        let left=xy[row.nodes[0]][0];let span=xy[*row.nodes.last().unwrap()][0]-left;
        for p in row.nodes.windows(2) {
            let mid=[0.5*(xy[p[0]][0]+xy[p[1]][0]),0.5*(xy[p[0]][1]+xy[p[1]][1])];
            let u=(mid[0]-left)/span;
            let h=max_h*0.001*(0.30+0.70*(4.0*u*(1.0-u)).min(1.0));
            beams.push(Beam{a:p[0],b:p[1],width:0.025,height:h,z:-0.5*(thickness(mid)+h),
                e:9.0e9,g:0.6e9,rho:400.0,group:format!("sugar_pine_rib_{}",17-index)});
        }
    }
    for k in 0..3 {
        for pair in rows.windows(2) {
            if let (Some(a),Some(b))=(pair[0].features[k],pair[1].features[k]) {
                let mid=[0.5*(xy[a][0]+xy[b][0]),0.5*(xy[a][1]+xy[b][1])];
                let s=0.5*(pair[0].station+pair[1].station);
                let first=station(curves[k][0]);let last=station(*curves[k].last().unwrap());
                let u=((s-first)/(last-first)).clamp(0.,1.);
                let (w,h,e,g,rho,name,side)=if k==2 {(0.040,0.035,9e9,0.6e9,400.,"cutoff_bar",-1.)}
                    else {(0.026,0.044-0.022*u,12.6e9,0.9e9,705.,if k==0{"maple_long_bridge"}else{"maple_bass_bridge"},1.)};
                beams.push(Beam{a,b,width:w,height:h,z:side*0.5*(thickness(mid)+h),e,g,rho,group:name.into()});
            }
        }
    }
    for beam in &beams {
        if !edges.contains_key(&(beam.a.min(beam.b),beam.a.max(beam.b))) {
            return Err("structural beam crosses a panel triangle interior".into());
        }
    }
    let contacts = if published_contacts {
        "RT-0425 Appendix A coupling points projected to reconstructed bridges; four extrapolated end keys"
    } else { "estimated per-key bridge stations" };
    let mut geometry=format!("frankensim-board-geometry-si-v1\nsource,mixed,Approximate Model D reconstruction: Boutillon Ege Paulello 2012 Fig 2/Table 1; Steinway Model D specifications; {contacts}; see steinway_d.rs for estimates\nsupport,clamped\npretension,0\ndamping,0.015\n");
    let mut obj=String::from("# Source-derived Model D SOUND BOARD assembly, SI metres, z up; not full piano CAD\n# Same flat structural image as the FSB. No invented crown or measured-data claim.\no sitka_spruce_panel\n");
    for (i,p) in xy.iter().enumerate() {
        writeln!(geometry,"node,{i},{:.17e},{:.17e}",p[0],p[1]).unwrap();
        writeln!(obj,"v {:.9} {:.9} {:.9}",p[0],p[1],0.5*thickness(*p)).unwrap();
    }
    for p in &xy {writeln!(obj,"v {:.9} {:.9} {:.9}",p[0],p[1],-0.5*thickness(*p)).unwrap();}
    let n=xy.len();
    // Grain direction in SI x/y, normal to a rib's (1,0.8) tangent.
    let grain=1.0f64.atan2(-SLOPE);
    for (i,&[a,b,c]) in tris.iter().enumerate() {
        let mid=[(xy[a][0]+xy[b][0]+xy[c][0])/3.,(xy[a][1]+xy[b][1]+xy[c][1])/3.];
        writeln!(geometry,"triangle,{i},{a},{b},{c},{:.17e},380,1.15e10,7.4e8,0.35,6.5e8,{grain:.17e}",thickness(mid)).unwrap();
        writeln!(obj,"f {} {} {}\nf {} {} {}",a+1,b+1,c+1,c+1+n,b+1+n,a+1+n).unwrap();
    }
    let mut fixed=std::collections::BTreeSet::new();
    for (&(a,b),&count) in &edges {if count==1 {fixed.insert(a);fixed.insert(b);}}
    for node in fixed {writeln!(geometry,"fixed,{node}").unwrap();}
    // OBJ side walls follow oriented boundary edges, keeping the panel closed.
    for &[a,b,c] in &tris {for (u,v) in [(a,b),(b,c),(c,a)] {
        if edges[&(u.min(v),u.max(v))]==1 {writeln!(obj,"f {} {} {} {}",v+1,u+1,u+1+n,v+1+n).unwrap();}
    }}
    let mut vertices=2*n;
    for beam in &beams {
        let (w,h)=(beam.width,beam.height);let area=w*h;let inertia=w*h.powi(3)/12.;
        // Saint-Venant torsion approximation for a rectangular section.
        let (long,short)=(w.max(h),w.min(h));let ratio=short/long;
        let torsion=long*short.powi(3)*(1./3.-0.21*ratio*(1.-ratio.powi(4)/12.));
        writeln!(geometry,"stiffener,{:.17e},{:.17e},{area:.17e},{inertia:.17e},{torsion:.17e},{:.17e},{:.17e},{},{}",beam.e,beam.g,beam.z,beam.rho,beam.a,beam.b).unwrap();
        writeln!(obj,"g {}",beam.group).unwrap();
        let (a,b)=(xy[beam.a],xy[beam.b]);let len=((b[0]-a[0]).powi(2)+(b[1]-a[1]).powi(2)).sqrt();
        let normal=[-(b[1]-a[1])/len*w/2.,(b[0]-a[0])/len*w/2.];
        for z in [beam.z-h/2.,beam.z+h/2.] {for (p,sign) in [(a,-1.),(b,-1.),(b,1.),(a,1.)] {
            writeln!(obj,"v {:.9} {:.9} {:.9}",p[0]+sign*normal[0],p[1]+sign*normal[1],z).unwrap();
        }}
        for face in [[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]] {
            writeln!(obj,"f {} {} {} {}",vertices+face[0]+1,vertices+face[1]+1,vertices+face[2]+1,vertices+face[3]+1).unwrap();
        }
        vertices+=8;
    }
    // The legacy preset estimated stations by arc length. The explicit
    // RT-0425 variant places the 84 published contacts on the closest point of
    // each reconstructed bridge. Four end keys remain extrapolations.
    for key in 21u8..=108 {
        let p=if published_contacts { published_bridge_point(key).0 }
            else {
                let raster=if key<=40 {curve_point(BASS_BRIDGE,0.95-0.90*f64::from(key-21)/19.)}
                    else {curve_point(LONG_BRIDGE,0.04+0.92*f64::from(key-41)/67.)};
                si(raster)
            };
        let (tri,w)=locate(p,&xy,&tris)?;
        writeln!(geometry,"bridge,{key},{tri},{:.17e},{:.17e},{:.17e}",w[0],w[1],w[2]).unwrap();
    }
    Ok(Preset{geometry,obj})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_reconstruction_is_admitted_by_existing_plate_front_door() {
        let p=build(6).unwrap();
        super::super::board_geometry::BoardGeometry::read(&p.geometry).unwrap();
        assert_eq!(p.geometry.lines().filter(|l|l.starts_with("bridge,")).count(),88);
        assert!(p.obj.contains("sugar_pine_rib_17"));
        assert!(p.obj.contains("sugar_pine_rib_1\n"));
        assert!(p.obj.contains("maple_bass_bridge"));
        assert!(p.obj.contains("cutoff_bar"));
        assert_eq!(p.geometry,build(6).unwrap().geometry);
    }
    #[test]
    fn generated_bridge_contacts_lie_on_conforming_structural_edges() {
        for divisions in [6, 12, 24, 28, 32] {
            for published in [false, true] {
                let p=build_with_contacts(divisions,published).unwrap();
                super::super::board_geometry::BoardGeometry::read(&p.geometry).unwrap();
                let mut edges=std::collections::BTreeSet::new();
                for row in p.geometry.lines().filter(|row|row.starts_with("triangle,")) {
                    let nodes=row.split(',').skip(2).take(3)
                        .map(|x|x.parse::<usize>().unwrap()).collect::<Vec<_>>();
                    for pair in [(nodes[0],nodes[1]),(nodes[1],nodes[2]),(nodes[2],nodes[0])] {
                        edges.insert((pair.0.min(pair.1),pair.0.max(pair.1)));
                    }
                }
                let mut beam_count=0;
                for row in p.geometry.lines().filter(|row|row.starts_with("stiffener,")) {
                    let nodes=row.split(',').rev().take(2)
                        .map(|x|x.parse::<usize>().unwrap()).collect::<Vec<_>>();
                    assert!(edges.contains(&(nodes[0].min(nodes[1]),nodes[0].max(nodes[1]))),
                        "mesh={divisions} published={published}: beam crosses panel: {row}");
                    beam_count+=1;
                }
                assert!(beam_count>0);
                let mut count=0;
                for row in p.geometry.lines().filter(|row|row.starts_with("bridge,")) {
                    let weights=row.split(',').skip(3).map(|x|x.parse::<f64>().unwrap()).collect::<Vec<_>>();
                    assert_eq!(weights.len(),3);
                    assert!(weights.iter().any(|w|w.abs()<1e-8),
                        "mesh={divisions} published={published}: contact inside triangle: {row}");
                    count+=1;
                }
                assert_eq!(count,88);
            }
        }
    }
    #[test]
    fn source_scale_and_bounds_are_physical() {
        let spacing=(LAST_RIB-FIRST_RIB)/16.*scale()/(1.+SLOPE*SLOPE).sqrt();
        assert!((spacing-SPACING_M).abs()<1e-14);
        for p in OUTLINE {let q=si(*p);assert!((0.0..1.56).contains(&q[0]));assert!((0.0..2.74).contains(&q[1]));}
        assert_eq!(RIBS.len(),17);
        assert!(build(0).is_err());assert!(build(33).is_err());
        let fine=build(32).unwrap();
        super::super::board_geometry::BoardGeometry::read(&fine.geometry).unwrap();
        assert_eq!(fine.geometry.lines().filter(|l|l.starts_with("bridge,")).count(),88);
    }
    #[test]
    fn rt0425_contacts_cover_published_keys_and_remain_on_reconstructed_bridges() {
        let p = build_with_rt0425_contacts(12).unwrap();
        super::super::board_geometry::BoardGeometry::read(&p.geometry).unwrap();
        assert_eq!(p.geometry.lines().filter(|l| l.starts_with("bridge,")).count(), 88);
        for key in 24..=107 {
            let (_, correction) = published_bridge_point(key);
            assert!(correction < 0.08, "key {key}: {correction} m projection");
        }
        let (a4, correction) = published_bridge_point(69);
        assert!((a4[0] - 0.81).abs() < 0.03);
        assert!((a4[1] - 0.30).abs() < 0.03);
        assert!(correction < 0.03);
        assert_ne!(p.geometry, build(12).unwrap().geometry);
    }
}
