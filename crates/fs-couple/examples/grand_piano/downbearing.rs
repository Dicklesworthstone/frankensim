//! Declared per-course static downbearing -> the existing full shell operator.
//! Loads are dead vertical forces at the SAME barycentric bridge/arm ports as
//! the dynamics. Positive input means downward N, total over the whole course.
//! The supplied crown must explicitly be an UNLOADED reference. A measured
//! already-loaded shape must not receive the same downbearing a second time.
use super::{CrownedBoard, Site, ShellMesh, ShellModel, ShellSupport,
    assemble_stiffened_shell, cross, meaningful, number, MAX_HEIGHT_M};
use fs_plate::shell::preload::{PreloadOptions, PreloadProblem, PreloadedShell,
    equilibrate_stiffened_shell};
use fs_solver::{NewtonKrylovConfig, NewtonKrylovState, NonlinearProblem};
use std::collections::BTreeMap;

pub(super) struct Specification {
    source: String,
    forces: BTreeMap<u8,f64>,
}
impl Specification {
    pub(super) fn read(text:&str)->Result<Option<Self>,String> {
        let mut unloaded=false;let mut source=None;let mut forces=BTreeMap::new();
        for row in meaningful(text) {
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            match f[0] {
                "preload-reference"=>{
                    if f.len()!=2 || f[1]!="unloaded" || unloaded {
                        return Err("declare preload-reference,unloaded exactly once; do not reload a measured loaded shape".into());
                    }
                    unloaded=true;
                }
                "downbearing-source"=>{
                    if f.len()<3 || !["estimated","mixed","published","measured"].contains(&f[1])
                        || f[2..].iter().all(|s|s.is_empty()) || row.len()>8192 || source.is_some() {
                        return Err("downbearing needs one bounded source authority and attribution".into());
                    }
                    source=Some(f[1..].join(","));
                }
                "downbearing"=>{
                    if f.len()!=3 {return Err("expected downbearing,key,total_downward_force_n".into());}
                    let key:u8=f[1].parse().map_err(|_|"invalid downbearing key")?;
                    let force=number(f[2])?;
                    if !(21..=108).contains(&key) || !(0.0..=10_000.0).contains(&force)
                        || forces.insert(key,force).is_some() {
                        return Err("downbearing key must be unique, force finite in [0,10000] N per course".into());
                    }
                }
                _=>{}
            }
        }
        if !unloaded && source.is_none() && forces.is_empty() {return Ok(None);}
        if !unloaded || source.is_none() || forces.is_empty() {
            return Err("downbearing requires unloaded reference, attributed source and explicit course loads".into());
        }
        Ok(Some(Self {source:source.expect("checked source"),forces}))
    }
    pub(super) fn check_sites(&self,sites:&[Site])->Result<(),String> {
        if self.forces.len()!=sites.len() || sites.iter().any(|s|!self.forces.contains_key(&s.key)) {
            return Err("downbearing must cover every geometric bridge station, including unplayed keys; use explicit zero rows for unloaded courses".into());
        }
        Ok(())
    }
    fn loads(&self,board:&CrownedBoard)->Result<Vec<f64>,String> {
        self.check_sites(&board.sites)?;
        let mut loads=vec![0.;6*board.mesh.nodes.len()];
        for site in &board.sites {
            let force=[0.,0.,-self.forces[&site.key]];
            let moment=cross(site.arm,force);
            for (a,&node) in board.mesh.tris[site.tri].iter().enumerate() {for c in 0..3 {
                loads[6*node+c]+=site.weights[a]*force[c];
                loads[6*node+3+c]+=site.weights[a]*moment[c];
            }}
        }
        if loads.iter().any(|v|!v.is_finite()) {return Err("bridge load projection overflow".into());}
        Ok(loads)
    }
}

// fs-solver owns Newton/Krylov/globalization. Its dependency graph includes
// FEEC and optionally fs-couple, so inject it from this example rather than
// introducing a runtime dependency cycle in the generic fs-plate owner.
struct Owner<'a,'b>(&'a PreloadProblem<'b>);
impl NonlinearProblem for Owner<'_,'_> {
    fn dimension(&self)->usize {self.0.dimension()}
    fn residual(&self,x:&[f64],out:&mut[f64]) {self.0.residual(x,out);}
    fn jacobian_apply(&self,x:&[f64],direction:&[f64],out:&mut[f64]) {
        self.0.jacobian_apply(x,direction,out);
    }
}
fn solve(problem:&PreloadProblem<'_>,x:Vec<f64>,limit:usize)->Result<(Vec<f64>,usize),String> {
    let owner=Owner(problem);
    let config=NewtonKrylovConfig {absolute_tolerance:1e-14,relative_tolerance:1e-11,
        forcing_maximum:0.1,..NewtonKrylovConfig::default()};
    let mut state=NewtonKrylovState::new(&owner,x,config).map_err(|e|e.to_string())?;
    let report=state.run(&owner,limit);
    if !report.converged {
        return Err(format!("downbearing Newton did not converge: {:?}, scaled residual {}",
            report.diagnosis,report.residual_norm));
    }
    Ok((state.x,report.iterations))
}
fn equilibrium(board:&CrownedBoard,spec:&Specification)->Result<PreloadedShell,String> {
    equilibrate_stiffened_shell(&board.mesh,&board.sections,&board.fixed,ShellSupport::Clamped,
        &board.beams,&spec.loads(board)?,PreloadOptions::default(),&mut solve)
}

/// Return the actual small-signal pencil and the loaded observation geometry.
/// M and linear beam arms retain reference values. The acoustic projection is
/// still a flat-baffle approximation, now using equilibrium normals and area.
pub(super) fn prepare(board:&CrownedBoard)->Result<(ShellModel,ShellMesh,String),String> {
    let Some(spec)=&board.preload else {
        let model=assemble_stiffened_shell(&board.mesh,&board.sections,&board.fixed,
            ShellSupport::Clamped,&board.beams).map_err(|e|e.to_string())?;
        return Ok((model,board.mesh.clone(),"no downbearing equilibrium".into()));
    };
    let report=equilibrium(board,spec)?;
    let mut mesh=board.mesh.clone();let mut maximum_displacement=0.0_f64;
    for (node,p) in mesh.nodes.iter_mut().enumerate() {
        let u=&report.displacement[6*node..6*node+3];
        maximum_displacement=maximum_displacement.max(u.iter().fold(0.0_f64,|n,v|n.hypot(*v)));
        for c in 0..3 {p[c]+=u[c];}
        if p.iter().any(|v|!v.is_finite()) || p[2].abs()>MAX_HEIGHT_M {
            return Err("loaded board leaves the admitted shallow-board height budget".into());
        }
    }
    for e in 0..mesh.tris.len() {
        if mesh.facet(e).map_err(|e|e.to_string())?.frame[2][2]<0.95 {
            return Err("loaded board is not an upward shallow acoustic graph".into());
        }
    }
    let summary=format!("nonlinear dead-load downbearing equilibrium: source [{}], total {} N, maximum translation {} m, physical residual {} N (moment norm length {} m), {} Newton iterations / {} evaluations, stored static energy {} J; reference mass and linear beams; no follower loads or beam buckling",
        spec.source,spec.forces.values().sum::<f64>(),maximum_displacement,report.residual_force_n,
        report.norm_length_m,report.iterations,report.evaluations,report.stored_energy_j);
    Ok((report.model,mesh,summary))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write;
    // Same small authored crown as the parent's structural regression, not
    // a measured Steinway; actual per-key source strings supply the last test.
    fn text(force:Option<f64>)->String {
        let mut f=String::from("frankensim-board-geometry-si-v1\nsource,estimated,downbearing regression\nsupport,clamped\npretension,0\ndamping,0.01\n");
        for (i,p) in [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.5,0.5]].iter().enumerate() {
            writeln!(f,"node,{i},{},{}",p[0],p[1]).unwrap();
        }
        for (i,t) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
            writeln!(f,"triangle,{i},{},{},{},0.008,450,1e10,8e8,0.3,6e8,0.27",t[0],t[1],t[2]).unwrap();
        }
        f.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
        let mut out=super::super::elevate(&f,&[0.,0.,0.,0.,0.015],"unloaded authored crown").unwrap();
        if let Some(force)=force {writeln!(out,"preload-reference,unloaded\ndownbearing-source,estimated,regression load\ndownbearing,69,{force}").unwrap();}
        out
    }
    #[test]
    fn downbearing_changes_the_actual_eigenproblem_and_keeps_static_force_balance() {
        let before=CrownedBoard::read(&text(None)).unwrap();
        let after=CrownedBoard::read(&text(Some(10.))).unwrap();
        let report=equilibrium(&after,after.preload.as_ref().unwrap()).unwrap();
        assert!(report.displacement[26]<0.);assert!(report.residual_force_n<1e-6);
        assert!(report.iterations>0);assert!(report.stored_energy_j>0.);
        let a=before.prepare(&[69],400.).unwrap();let b=after.prepare(&[69],400.).unwrap();
        assert!((a.modes[0].frequency_hz-b.modes[0].frequency_hz).abs()>0.01);
        assert_eq!(a.mass_kg,b.mass_kg);
        assert_ne!(a.surface[0].area_m2,b.surface[0].area_m2);
        assert!(b.provenance.contains("physical residual"));
    }
    #[test]
    fn bridge_load_and_motion_are_work_conjugate_with_a_lever_arm() {
        let mut board=CrownedBoard::read(&text(Some(10.))).unwrap();
        board.sites[0].weights=[0.2,0.3,0.5];board.sites[0].arm=[0.02,-0.03,0.04];
        let load=board.preload.as_ref().unwrap().loads(&board).unwrap();
        let trial:Vec<f64>=(0..30).map(|i|(i as f64*0.71).sin()*0.001).collect();
        let work:f64=trial.iter().zip(&load).map(|(a,b)|a*b).sum();
        let site=&board.sites[0];let mut dz=0.;
        for (i,&node) in board.mesh.tris[site.tri].iter().enumerate() {
            let u=[trial[6*node],trial[6*node+1],trial[6*node+2]];
            let theta=[trial[6*node+3],trial[6*node+4],trial[6*node+5]];
            dz+=site.weights[i]*(u[2]+cross(theta,site.arm)[2]);
        }
        assert!((work+10.*dz).abs()<1e-16);
        assert!((load.chunks_exact(6).map(|p|p[2]).sum::<f64>()+10.).abs()<1e-14);
    }
    #[test]
    fn loading_requires_unloaded_geometry_and_every_station_not_only_played_keys() {
        let good=text(Some(10.));
        for bad in [good.replace("preload-reference,unloaded\n",""),
            good.replace("preload-reference,unloaded","preload-reference,loaded"),
            good.replace("downbearing-source,estimated,regression load\n",""),
            good.replace("downbearing,69,10","downbearing,69,NaN"),
            good.replace("downbearing,69,10","downbearing,69,-1"),
            format!("{good}downbearing,69,10\n"),format!("{good}bridge,60,1,0,0,1\n"),
            format!("{good}downbearing,60,0\n")] {assert!(CrownedBoard::read(&bad).is_err());}
        let complete=format!("{good}bridge,60,1,0,0,1\ndownbearing,60,5\n");
        let b=CrownedBoard::read(&complete).unwrap();
        let load=b.preload.as_ref().unwrap().loads(&b).unwrap();
        assert_eq!(load.chunks_exact(6).map(|p|p[2]).sum::<f64>(),-15.);
    }
    #[test]
    fn zero_downbearing_preserves_the_native_pencil_and_pressure_surface() {
        let a=CrownedBoard::read(&text(None)).unwrap().prepare(&[69],400.).unwrap();
        let b=CrownedBoard::read(&text(Some(0.))).unwrap().prepare(&[69],400.).unwrap();
        assert_eq!(a.modes.len(),b.modes.len());
        for (a,b) in a.modes.iter().zip(&b.modes) {
            assert_eq!(a.frequency_hz,b.frequency_hz);assert_eq!(a.bridge,b.bridge);assert_eq!(a.volume,b.volume);
        }
        for (a,b) in a.surface.iter().zip(&b.surface) {
            assert_eq!(a.position_m,b.position_m);assert_eq!(a.area_m2,b.area_m2);assert_eq!(a.mode_shape,b.mode_shape);
        }
    }
    #[test]
    fn source_hammer_mechanics_reach_pressure_through_the_loaded_board() {
        use super::super::super::{steinway_scale,engine,performance,audio};
        let course=steinway_scale::courses().unwrap().into_iter().find(|c|c.midi==69).unwrap();
        let render=|force| {
            let board=CrownedBoard::read(&text(force)).unwrap().prepare(&[69],400.).unwrap();
            let piano=engine::Instrument::new_with_course_shanks(vec![course],&board.modes,
                48_000,4,12,true,vec![steinway_scale::hammer_material(&course).unwrap()],engine::ShankGeometry::published()).unwrap();
            let score=performance::Performance::read("sample,event,key,value\n0,note_on,69,2\n600,note_off,69,0\n",&[69],1800).unwrap();
            let mut stream=audio::AudioStream::new(piano,score,Some(&board.surface),[0.5,0.5,1.],
                fs_bem::helmholtz::Medium::air(),1000.).unwrap();
            let mut output=vec![0.;1800];stream.render_block(&mut output).unwrap();
            assert!(output.iter().any(|v|v.abs()>1e-10));
            let instrument=stream.instrument();
            let residual=instrument.accounting.input_work_j-instrument.energy_j()-instrument.accounting.dissipated_j();
            assert!(residual.abs()<1e-7);
            output
        };
        assert_ne!(render(None),render(Some(10.)));
    }
}
