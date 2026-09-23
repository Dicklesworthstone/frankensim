use super::*;
use super::super::Specimen;
use super::super::super::{Stroke, splash_with_sticks, muffling};
use std::fmt::Write;
use std::path::Path;

fn reference_text() -> String {
    let spec = Specimen::reference(); let shell = spec.build().unwrap();
    let mut s = format!("{HEADER}\nband_hz,50,1200\nstrike,0.06,0.01\n");
    writeln!(s,"material,23,{:.17e},{:.17e},{:.17e}",spec.young_pa,spec.poisson,spec.density_kg_m3).unwrap();
    for (i,(p,h)) in shell.mesh.nodes.iter().zip(&shell.nodal_thickness_m).enumerate() {
        writeln!(s,"node,{},{:.17e},{:.17e},{:.17e},{:.17e}",100+3*i,p[0],p[1],p[2],h).unwrap();
    }
    for (i,t) in shell.mesh.tris.iter().enumerate() {
        writeln!(s,"triangle,{},{},{},{},23",10+7*i,100+3*t[0],100+3*t[1],100+3*t[2]).unwrap();
    }
    s
}

#[test]
fn sparse_reordered_records_preserve_every_shape_thickness_and_material() {
    let text = reference_text(); let (a,band,strike) = parse(&text).unwrap();
    let mut lines:Vec<_> = text.lines().skip(1).collect(); lines.reverse();
    let (b,b_band,b_strike) = parse(&format!("{HEADER}\n{}",lines.join("\n"))).unwrap();
    assert_eq!(a.mesh,b.mesh); assert_eq!(a.nodal_thickness_m,b.nodal_thickness_m);
    assert_eq!(a.mass_kg,b.mass_kg); assert_eq!(band,b_band); assert_eq!(strike,b_strike);
    let reference = Specimen::reference().build().unwrap();
    assert_eq!(a.mesh,reference.mesh); assert_eq!(a.nodal_thickness_m,reference.nodal_thickness_m);
    assert_eq!(a.mass_kg,reference.mass_kg);
    for (x,y) in a.sections.iter().zip(&reference.sections) { assert_eq!(x.d,y.d); }
    let input = Specimen::parse_mesh(&text).unwrap();
    assert_eq!(input.default_strike_position().unwrap(),[0.06,0.01]);
    assert_eq!(input.input_label(),"supplied 3D mesh");
}

#[test]
fn supplied_geometry_changes_mechanics_without_hidden_profile_or_material_fallback() {
    let text = reference_text(); let (old,_,_) = parse(&text).unwrap();
    let mut lines: Vec<_> = text.lines().map(String::from).collect();
    let node = lines.iter_mut().find(|l|l.starts_with("node,400,")).unwrap();
    let mut fields: Vec<_> = node.split(',').map(String::from).collect();
    fields[4] = format!("{:.17e}",fields[4].parse::<f64>().unwrap()+0.0001);
    fields[5] = format!("{:.17e}",fields[5].parse::<f64>().unwrap()*1.2); *node = fields.join(",");
    lines.push("material,99,90e9,0.3,7800".into());
    let face = lines.iter_mut().find(|l|l.starts_with("triangle,10,")).unwrap();
    *face = face.strip_suffix(",23").unwrap().to_string()+",99";
    let (changed,_,_) = parse(&lines.join("\n")).unwrap();
    assert_ne!(changed.mesh.nodes,old.mesh.nodes); assert_eq!(changed.mesh.tris,old.mesh.tris);
    assert_ne!(changed.mass_kg,old.mass_kg); assert_ne!(changed.sections[0].d,old.sections[0].d);
    let a = old.assemble(&[],fs_plate::shell::ShellSupport::Free).unwrap();
    let b = changed.assemble(&[],fs_plate::shell::ShellSupport::Free).unwrap();
    assert_ne!(a.k.get(0,0),b.k.get(0,0));
}

#[test]
fn explicit_mesh_drives_the_same_shared_two_stick_shell_and_radiation_as_its_profile() {
    let text = reference_text(); let stroke = Stroke { position_m:Some([0.06,0.01]),speed_m_s:0.8 };
    let second = Some(Stroke { position_m:Some([-0.05,0.02]),speed_m_s:0.8 });
    let m = [muffling::Muffler { surface:muffling::Surface::Shell,position_m:[0.075,0.],resistance_n_s_m:0.2 }];
    let mut a = splash_with_sticks(192,2e-6,true,stroke,Some(Specimen::reference()),&m,second).unwrap();
    let mut b = splash_with_sticks(192,2e-6,true,stroke,Some(Specimen::parse_mesh(&text).unwrap()),&m,second).unwrap();
    assert_eq!(a.force,b.force); assert_eq!(a.observer_a,b.observer_a);
    assert_eq!(a.second_stick.unwrap().coordinate,b.second_stick.unwrap().coordinate);
    assert_eq!(a.acoustics.as_ref().unwrap().state_modes(),b.acoustics.as_ref().unwrap().state_modes());
    let bounds = fs_couple::render::plate::impact::ImpactSubstepConfig { max_depth:8,max_attempts:511 };
    a.system = a.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds).unwrap();
    b.system = b.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds).unwrap();
    let gate = fs_exec::CancelGate::new_clock_free(); let initial = a.system.state().to_vec();
    for _ in 0..192 {
        let x=a.system.step(&a.force,&gate).unwrap(); let y=b.system.step(&b.force,&gate).unwrap();
        assert_eq!(a.system.state(),b.system.state()); assert_eq!(x.stored_energy_j,y.stored_energy_j);
        assert!(x.balance_residual_j.abs()<1e-7);
    }
    assert_ne!(a.system.state(),initial);
    assert!(a.system.state()[2..].iter().any(|v|*v!=0.0));
}

#[test]
fn malformed_geometry_missing_fields_wrong_front_doors_and_hole_strikes_refuse() {
    let text = reference_text();
    for bad in [text.replace("strike,0.06,0.01","strike,0,0"),
        text.replace("band_hz,50,1200","band_hz,1200,50"),text.replace("strike,0.06,0.01\n",""),
        text.replace("triangle,10,100,196,199,23","triangle,10,100,196,199,999"),
        format!("{text}material,23,1e9,0.3,8000\n"),format!("{text}material,99,1e9,0.3,8000\n"),
        format!("{text}node,100,0,0,0,0.001\n"),format!("{text}unknown,1\n")] {
        assert!(parse(&bad).is_err(),"{bad}");
    }
    let mut missing = text.lines().filter(|l|!l.starts_with("node,100,")).collect::<Vec<_>>().join("\n");
    assert!(parse(&missing).is_err());
    missing = text.lines().filter(|l|!l.starts_with("material,")).collect::<Vec<_>>().join("\n");
    assert!(parse(&missing).is_err());
    assert!(parse(&"x".repeat(MAX_BYTES+1)).is_err());
    let p=Path::new("does-not-exist");
    assert!(super::super::admit_selection(Some(p),Some(p),"splash").is_err());
    assert!(super::super::load_selection(Some(p),Some(p)).is_err());
    for cmd in ["drum","snare","drum-modal","unknown"] {
        assert!(super::super::admit_selection(None,Some(p),cmd).is_err());
    }
    for cmd in ["splash","splash-wav","splash-mic"] {
        super::super::admit_selection(None,Some(p),cmd).unwrap();
    }
    for bad in [vec!["--shell-mesh".into()],vec!["--shell-mesh".into(),"--analytic-newton".into()],
        vec!["--shell-mesh".into(),"a".into(),"--shell-mesh".into(),"b".into()]] {
        let mut args = bad.clone(); assert!(super::super::mesh_option(&mut args).is_err()); assert_eq!(args,bad);
    }
    let mut args=["splash-mic","--shell-mesh","scan.fss","--analytic-newton"].map(String::from).to_vec();
    assert_eq!(super::super::mesh_option(&mut args).unwrap().unwrap(),Path::new("scan.fss"));
    assert_eq!(args,["splash-mic","--analytic-newton"]);
}

#[test]
fn profile_export_roundtrips_all_physical_fields_without_running_an_eigensolve() {
    let spec = Specimen::reference(); let before = spec.build().unwrap();
    let text = spec.mesh_text().unwrap(); let imported = Specimen::parse_mesh(&text).unwrap();
    let after = imported.build().unwrap();
    assert_eq!(before.mesh,after.mesh); assert_eq!(before.mass_kg,after.mass_kg);
    assert_eq!(before.nodal_thickness_m,after.nodal_thickness_m);
    for (a,b) in before.sections.iter().zip(&after.sections) { assert_eq!(a.d,b.d); }
    assert_eq!(spec.band_hz,imported.band_hz);
    assert_eq!(spec.default_strike_position().unwrap(),imported.default_strike_position().unwrap());
    assert!(imported.mesh_text().is_err(),"mesh input must not be reinterpreted as a profile");
    assert!(!super::super::export_command(&["splash".into()]).unwrap());
    assert!(super::super::export_command(&["export-shell-mesh".into()]).is_err());
    assert!(super::super::export_command(&["export-shell-mesh".into(),"out".into(),"profile".into(),"extra".into()]).is_err());
}
