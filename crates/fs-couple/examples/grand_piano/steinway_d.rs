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
    if !(4..=24).contains(&divisions) { return Err("Model D mesh divisions must be 4..24".into()); }
    let curves=[LONG_BRIDGE,BASS_BRIDGE,CUTOFF];
    let mut levels:Vec<f64>=OUTLINE.iter().copied().map(station).collect();
    levels.extend(RIBS.iter().map(|r|r.0));
    for curve in curves { levels.extend(curve.iter().copied().map(station)); }
    levels.sort_by(f64::total_cmp); levels.dedup_by(|a,b|(*a-*b).abs()<1e-8);
    let mut extra=Vec::new();
    for p in levels.windows(2) { if p[1]-p[0]>55.0 {extra.push((p[0]+p[1])*0.5);} }
    levels.extend(extra); levels.sort_by(f64::total_cmp);
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
    // Merge the normalized abscissae of adjacent rows: a conforming strip
    // triangulation, including one-node end caps, with no T-junctions.
    for pair in rows.windows(2) {
        let (a,b)=(&pair[0].nodes,&pair[1].nodes);let(mut i,mut j)=(0,0);
        while i+1<a.len()||j+1<b.len() {
            let next=|ids:&[usize],k:usize| {
                if k+1==ids.len(){f64::INFINITY}else{
                    (xy[ids[k+1]][0]-xy[ids[0]][0])/(xy[*ids.last().unwrap()][0]-xy[ids[0]][0])
                }
            };
            if next(a,i)<=next(b,j) {append_triangle(&mut tris,&xy,a[i],a[i+1],b[j]);i+=1;}
            else {append_triangle(&mut tris,&xy,a[i],b[j+1],b[j]);j+=1;}
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
    let mut geometry=String::from("frankensim-board-geometry-si-v1\nsource,mixed,Approximate Model D reconstruction: Boutillon Ege Paulello 2012 Fig 2/Table 1; Steinway Model D specifications; see steinway_d.rs for estimates\nsupport,clamped\npretension,0\ndamping,0.015\n");
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
    // Individual key/string bridge stations are not labelled in the source.
    // This estimated assignment is separate from the diagram-derived curves.
    // Bass low -> high approaches the junction; long bridge low -> high proceeds
    // toward the keyboard. End margins avoid placing strings directly on the rim.
    for key in 21u8..=108 {
        let p=if key<=40 {curve_point(BASS_BRIDGE,0.95-0.90*f64::from(key-21)/19.)}
            else {curve_point(LONG_BRIDGE,0.04+0.92*f64::from(key-41)/67.)};
        let (tri,w)=locate(si(p),&xy,&tris)?;
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
    fn source_scale_and_bounds_are_physical() {
        let spacing=(LAST_RIB-FIRST_RIB)/16.*scale()/(1.+SLOPE*SLOPE).sqrt();
        assert!((spacing-SPACING_M).abs()<1e-14);
        for p in OUTLINE {let q=si(*p);assert!((0.0..1.56).contains(&q[0]));assert!((0.0..2.74).contains(&q[1]));}
        assert_eq!(RIBS.len(),17);
        assert!(build(0).is_err());assert!(build(25).is_err());
    }
}
