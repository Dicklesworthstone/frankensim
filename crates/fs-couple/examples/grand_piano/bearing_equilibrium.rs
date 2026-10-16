//! Tuning equilibrium at supplied FINAL member tensions, with directions
//! evaluated on the deflected bridge. This is not a fixed-rest-length string
//! transient: the tuner maintains the target scale tensions while the board
//! settles. It preserves the exact per-member tensions used by playback.
//!
//! Export converged force rows on the ORIGINAL unloaded shell. Re-solving those
//! frozen loads recovers the same loaded surface and its BARE structural tangent;
//! the piano bank, not this file, supplies dynamic string stiffness and mass.
use super::{CrownedBoard, geometry, meaningful, number};
use super::super::downbearing;
use fs_plate::shell::preload::{PreloadOptions, equilibrate_tethered_shell};
use fs_plate::shell::preload::tethers::ShellTether;
use fs_plate::ShellSupport;
use std::{collections::BTreeMap,fmt::Write};

/// A native crowned board that can be played through the ordinary piano path.
pub struct SettledBoard {
    /// Original unloaded geometry and newly converged per-course force rows.
    pub geometry: String,
    /// Full reference-to-equilibrium coordinates: translation [m], rotation [rad].
    pub displacement: Vec<f64>,
    /// Physical shell-versus-current-string force residual [N].
    pub residual_force_n: f64,
    /// Maximum nodal translation norm [m].
    pub maximum_translation_m: f64,
    /// Maximum absolute change in an outgoing string span [m].
    pub maximum_length_change_m: f64,
    /// Number of outgoing spans; two for every admitted string member.
    pub spans: usize,
}

/// Solve direction-consistent downbearing at the supplied final tuned tensions.
/// Input admission and authority are shared with `bearing_geometry::apply`.
/// The old dead-reference path remains available and is not silently changed.
/// Refuse when the frozen-load playback preparation cannot recover the same
/// shallow loaded surface, including an unsupported/unstable bare board.
/// No input coordinates, scale cards, anchor points or material data mutate.
pub fn settle(board_text:&str,scale_text:&str,points_text:&str)->Result<SettledBoard,String> {
    // This front door already owns complete member/key coverage, SI/byte
    // limits, duplicate/unknown rows, source attribution and anchor distances.
    super::apply(board_text,scale_text,points_text)?;
    let board=CrownedBoard::read(board_text)?;
    let courses=geometry::read_scale(scale_text)?;
    let mut connections=BTreeMap::new();let mut source=None;
    for row in meaningful(points_text) {
        let f:Vec<_>=row.split(',').map(str::trim).collect();
        if f[0]=="source" {source=Some(f[1..].join(","));}
        if f[0]!="bearing" {continue;}
        let key:u8=f[1].parse().map_err(|_|"invalid bearing key")?;
        let member:usize=f[2].parse().map_err(|_|"invalid bearing member")?;
        let course=courses.iter().find(|c|c.midi==key).ok_or("bearing key absent from scale")?;
        let site=board.sites.iter().find(|s|s.key==key).ok_or("bearing key absent from board")?;
        let cents=(member as f64-0.5*(course.unison-1) as f64)*course.detune_cents;
        let tension=course.tension_at_cents(cents)?;
        for (side,start) in [3,6].into_iter().enumerate() {
            connections.insert((key,member,side),ShellTether {
                nodes:board.mesh.tris[site.tri],weights:site.weights,arm_m:site.arm,
                anchor_m:[number(f[start])?,number(f[start+1])?,number(f[start+2])?],
                reference_tension_n:tension,axial_stiffness_n_per_m:0.,
            });
        }
    }
    let keys:Vec<_>=connections.keys().map(|(key,_,_)|*key).collect();
    let connections:Vec<_>=connections.into_values().collect();
    let report=equilibrate_tethered_shell(&board.mesh,&board.sections,&board.fixed,
        ShellSupport::Clamped,&board.beams,&vec![0.;6*board.mesh.nodes.len()],&connections,
        PreloadOptions::default(),&mut downbearing::solve)?;
    let mut forces=BTreeMap::<u8,[f64;3]>::new();
    for (&key,response) in keys.iter().zip(&report.tether_responses) {
        let total=forces.entry(key).or_insert([0.;3]);
        for c in 0..3 {total[c]+=response.force_n[c];}
    }
    let mut out=format!("{}\npreload-reference,unloaded\ndownbearing-source,mixed,Converged bridge directions at prescribed final scale member tensions; support source [{}]; export retains bare-shell tangent for dynamic string coupling\n",
        board_text.trim_end(),source.ok_or("missing bearing attribution")?);
    for (key,force) in forces {
        writeln!(out,"bridge-load,{key},{:.17e},{:.17e},{:.17e}",force[0],force[1],force[2]).unwrap();
    }
    // Freeze ONLY the converged forces. Do not bake the tether Hessian into
    // the board and then let Bank add those same strings for a second time.
    // Use the real playback preparation to check the stable bare equilibrium
    // and acoustic geometry, not just serialization of plausible force rows.
    let frozen=CrownedBoard::read(&out)?;
    let (_,recovered,_)=downbearing::prepare(&frozen)?;
    let maximum_translation_m=report.displacement.chunks_exact(6).map(|u|
        u[..3].iter().fold(0.0_f64,|n,v|n.hypot(*v))).fold(0.0_f64,f64::max);
    let tolerance=1e-9+1e-6*maximum_translation_m;
    for (node,(a,b)) in board.mesh.nodes.iter().zip(&recovered.nodes).enumerate() {
        let error=(0..3).map(|c|(b[c]-a[c]-report.displacement[6*node+c]).abs()).fold(0.0_f64,f64::max);
        if !error.is_finite() || error>tolerance {
            return Err("frozen final forces do not recover the coupled tuning equilibrium".into());
        }
    }
    let maximum_length_change_m=report.tether_responses.iter()
        .map(|s|s.length_change_m.abs()).fold(0.0_f64,f64::max);
    Ok(SettledBoard {geometry:out,displacement:report.displacement,
        residual_force_n:report.residual_force_n,maximum_translation_m,
        maximum_length_change_m,spans:connections.len()})
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::input;
    #[test]
    fn final_bridge_forces_follow_the_solved_geometry_not_the_unloaded_chord() {
        let (board,scale,points)=input(2);
        let initial=super::super::apply(&board,&scale,&points).unwrap();
        let result=settle(&board,&scale,&points).unwrap();
        assert_eq!(result.spans,4);assert!(result.residual_force_n<1e-5);
        assert!(result.maximum_translation_m>1e-6);assert!(result.maximum_length_change_m>1e-8);
        let frozen=CrownedBoard::read(&result.geometry).unwrap();
        let original=CrownedBoard::read(&board).unwrap();
        assert_eq!(frozen.mesh,original.mesh,"do not erase prestress by changing reference shape");
        let loads=frozen.preload.as_ref().unwrap().loads(&frozen).unwrap();
        let p:[f64;3]=std::array::from_fn(|c|original.mesh.nodes[4][c]+result.displacement[24+c]);
        let mut expected=[0.;3];
        for anchor in [[0.,0.5,0.],[1.,0.5,0.]] {
            let v:[f64;3]=std::array::from_fn(|c|anchor[c]-p[c]);
            let length=v.iter().fold(0.0_f64,|n,v|n.hypot(*v));
            for c in 0..3 {expected[c]+=1400.*v[c]/length;}
        }
        for c in 0..3 {assert!((loads[24+c]-expected[c]).abs()<1e-8);}
        let before=CrownedBoard::read(&initial).unwrap();
        let before_load=before.preload.as_ref().unwrap().loads(&before).unwrap();
        assert!((before_load[26]-loads[26]).abs()>1e-4);
        assert_ne!(before.prepare(&[69],400.).unwrap().modes[0].frequency_hz,
            frozen.prepare(&[69],400.).unwrap().modes[0].frequency_hz);
    }
    #[test]
    fn balanced_tensions_preserve_zero_equilibrium_and_existing_admission() {
        let (board,scale,points)=input(1);
        let horizontal=points.replace("0,0.5,0,1,0.5,0","0,0.5,0.01,1,0.5,0.01");
        let result=settle(&board,&scale,&horizontal).unwrap();
        assert_eq!(result.maximum_translation_m,0.);assert_eq!(result.maximum_length_change_m,0.);
        assert!(settle(&result.geometry,&scale,&horizontal).is_err());
        assert!(settle(&board,&scale,&horizontal.replace("reference,unloaded\n","")).is_err());
    }
    #[test]
    fn settled_scale_tensions_and_loaded_modes_reach_the_existing_pressure_renderer() {
        use super::super::super::super::{engine,steinway_scale,performance,audio};
        let (board,scale,points)=input(1);
        let result=settle(&board,&scale,&points).unwrap();
        let initial=super::super::apply(&board,&scale,&points).unwrap();
        let course=geometry::read_scale(&scale).unwrap()[0];
        let render=|text:&str| {
            let b=CrownedBoard::read(text).unwrap().prepare(&[69],400.).unwrap();
            let piano=engine::Instrument::new_with_course_shanks(vec![course],&b.modes,
                48_000,4,12,true,vec![steinway_scale::hammer_material(&course).unwrap()],
                engine::ShankGeometry::published()).unwrap();
            let score=performance::Performance::read("sample,event,key,value\n0,note_on,69,2\n600,note_off,69,0\n",&[69],1800).unwrap();
            let mut stream=audio::AudioStream::new(piano,score,Some(&b.surface),[0.5,0.5,1.],
                fs_bem::helmholtz::Medium::air(),1000.).unwrap();
            let mut out=vec![0.;1800];stream.render_block(&mut out).unwrap();
            assert!(out.iter().any(|v|v.abs()>1e-10));
            let p=stream.instrument();
            assert!((p.accounting.input_work_j-p.energy_j()-p.accounting.dissipated_j()).abs()<1e-7);
            out
        };
        assert_ne!(render(&initial),render(&result.geometry));
    }
}
