//! Supplied string support coordinates + existing tension cards -> dead loads.
//! This converts geometry, not a note number or a desired frequency, into force.
//! Force is evaluated at the declared UNLOADED reference and then held during
//! equilibrium. No follower-load/axial-stretch update or bridge-pin friction.
use super::{CrownedBoard, meaningful, number, MAX_BYTES};
use super::super::geometry;
use std::{collections::{BTreeMap,BTreeSet},fmt::Write};

pub const HEADER:&str="frankensim-bearing-points-si-v1";
pub const MAX_BEARING_BYTES:usize=512*1024;

/// `bearing,key,member,front_x,front_y,front_z,rear_x,rear_y,rear_z` in SI metres
/// in the SAME Cartesian chart as the crowned board. Member is zero-based and
/// every string of every course must be supplied. The two points are the next
/// fixed/contact supports on either side of that course's bridge bearing, not
/// necessarily tuning/hitch pins. The bridge position is barycentric + arm.
/// Material/source metadata remain explicit; no measured provenance is invented.
pub fn apply(board_text:&str,scale_text:&str,points_text:&str)->Result<String,String> {
    if board_text.len()>MAX_BYTES || scale_text.len()>MAX_BEARING_BYTES
        || points_text.len()>MAX_BEARING_BYTES {return Err("bearing/scale input exceeds cold byte budget".into());}
    let board=CrownedBoard::read(board_text)?;
    if board.preload.is_some() {return Err("input already has loads; use its original unloaded reference".into());}
    let courses=geometry::read_scale(scale_text)?;
    if courses.len()!=board.sites.len() || courses.iter().any(|c|!board.sites.iter().any(|s|s.key==c.midi)) {
        return Err("string scale must cover exactly the geometric bridge keys".into());
    }
    let mut rows=meaningful(points_text);
    if rows.next()!=Some(HEADER) {return Err(format!("expected {HEADER}"));}
    let mut source=None;let mut unloaded=false;let mut seen=BTreeSet::new();
    let mut forces=BTreeMap::<u8,[f64;3]>::new();
    for row in rows {
        let f:Vec<_>=row.split(',').map(str::trim).collect();
        match f[0] {
            "source" if f.len()>=3=>{
                if source.is_some() || row.len()>8192
                    || !["measured","published","estimated","mixed"].contains(&f[1])
                    || f[2..].iter().all(|x|x.is_empty()) {
                    return Err("bearing geometry requires one attributed source row".into());
                }
                source=Some(f[1..].join(","));
            }
            "reference" if f.len()==2 && f[1]=="unloaded" && !unloaded=>unloaded=true,
            "bearing" if f.len()==9=>{
                let key:u8=f[1].parse().map_err(|_|"invalid bearing key")?;
                let member:usize=f[2].parse().map_err(|_|"invalid bearing member")?;
                let course=courses.iter().find(|c|c.midi==key).ok_or("bearing key absent from scale")?;
                if member>=course.unison || !seen.insert((key,member)) {
                    return Err("missing/out-of-range or duplicate unison bearing member".into());
                }
                let site=board.sites.iter().find(|s|s.key==key).ok_or("bearing key absent from geometry")?;
                let point:[f64;3]=std::array::from_fn(|c| site.arm[c]+(0..3).map(|i|
                    site.weights[i]*board.mesh.nodes[board.mesh.tris[site.tri][i]][c]).sum::<f64>());
                let cents=(member as f64-0.5*(course.unison-1) as f64)*course.detune_cents;
                let tension=course.tension_at_cents(cents)?;
                let force=forces.entry(key).or_insert([0.;3]);
                for start in [3,6] {
                    let support=[number(f[start])?,number(f[start+1])?,number(f[start+2])?];
                    let direction:[f64;3]=std::array::from_fn(|c|support[c]-point[c]);
                    let length=direction.iter().fold(0.0_f64,|n,v|n.hypot(*v));
                    if !length.is_finite() || !(1e-5..=10.0).contains(&length) {
                        return Err("bearing support must be 10 micrometres to 10 metres from the bridge point".into());
                    }
                    for c in 0..3 {force[c]+=tension*direction[c]/length;}
                }
            }
            _=>return Err("unknown bearing row or field count; reference must be unloaded".into()),
        }
    }
    if !unloaded || source.is_none() || seen.len()!=courses.iter().map(|c|c.unison).sum::<usize>() {
        return Err("bearing file must declare unloaded reference, attributed source, and both supports of every string".into());
    }
    // The union is mixed authority: tension cards and reference geometry need
    // not share the bearing source's epistemic status. Retain its attribution.
    let mut out=format!("{}\npreload-reference,unloaded\ndownbearing-source,mixed,Forces from supplied scale member tensions and support coordinates [{}]; reference-frame dead forces, no follower update\n",board_text.trim_end(),source.unwrap());
    for (key,force) in forces {
        writeln!(out,"bridge-load,{key},{:.17e},{:.17e},{:.17e}",force[0],force[1],force[2]).unwrap();
    }
    CrownedBoard::read(&out)?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn input(unison:usize)->(String,String,String) {
        let flat="frankensim-board-geometry-si-v1\nsource,estimated,anchor regression\nsupport,clamped\npretension,0\ndamping,0.01\nnode,0,0,0\nnode,1,1,0\nnode,2,1,1\nnode,3,0,1\nnode,4,0.5,0.5\ntriangle,0,0,1,4,0.008,450,1e10,8e8,0.3,6e8,0.27\ntriangle,1,1,2,4,0.008,450,1e10,8e8,0.3,6e8,0.27\ntriangle,2,2,3,4,0.008,450,1e10,8e8,0.3,6e8,0.27\ntriangle,3,3,0,4,0.008,450,1e10,8e8,0.3,6e8,0.27\nfixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n";
        let board=super::super::elevate(flat,&[0.,0.,0.,0.,0.01],"authored reference").unwrap();
        let mut course=geometry::demonstration_scale().unwrap()[48];
        course.unison=unison;course.tension_n=700.;course.detune_cents=0.;
        let scale=geometry::write_scale(&[course]);
        let mut points=format!("{HEADER}\nsource,estimated,3-D point regression\nreference,unloaded\n");
        for member in 0..unison {writeln!(points,"bearing,69,{member},0,0.5,0,1,0.5,0").unwrap();}
        (board,scale,points)
    }
    #[test]
    fn support_geometry_and_member_tensions_produce_the_vector_resultant() {
        for members in [1,2,3] {
            let (board,scale,points)=input(members);
            let result=apply(&board,&scale,&points).unwrap();
            let b=CrownedBoard::read(&result).unwrap();
            let loads=b.preload.as_ref().unwrap().loads(&b).unwrap();
            let force: [f64;3]=std::array::from_fn(|c|loads.chunks_exact(6).map(|p|p[c]).sum());
            let expected=-(members as f64)*2.*700.*0.01/(0.5_f64.powi(2)+0.01_f64.powi(2)).sqrt();
            assert!(force[0].abs()<1e-10);assert_eq!(force[1],0.);
            assert!((force[2]-expected).abs()<1e-10);
        }
    }
    #[test]
    fn asymmetric_supports_create_in_plane_force_and_actual_bridge_torque() {
        let (mut board,scale,mut points)=input(1);
        board.push_str("bridge_arm,69,0,0,0.02\n");
        points=points.replace("1,0.5,0","1,0.55,0");
        let result=apply(&board,&scale,&points).unwrap();let b=CrownedBoard::read(&result).unwrap();
        let loads=b.preload.as_ref().unwrap().loads(&b).unwrap();
        let center=&loads[24..30];
        assert!(center[1]>0.);assert!(center[2]<0.);
        assert!((center[3]+0.02*center[1]).abs()<1e-12);
        assert!((center[4]-0.02*center[0]).abs()<1e-12);
    }
    #[test]
    fn missing_members_and_already_loaded_geometry_never_acquire_default_anchors() {
        let (board,scale,points)=input(2);
        let result=apply(&board,&scale,&points).unwrap();
        assert!(apply(&result,&scale,&points).is_err());
        for bad in [points.replace("reference,unloaded\n",""),
            points.replace("bearing,69,1,0,0.5,0,1,0.5,0\n",""),
            points.replace("bearing,69,1","bearing,69,0"),
            points.replace("bearing,69,1","bearing,69,2"),
            points.replace("0,0.5,0,1,0.5,0","0.5,0.5,0.01,1,0.5,0"),
            points.replace("bearing,69,1","bearing,60,1")] {
            assert!(apply(&board,&scale,&bad).is_err());
        }
    }
}
