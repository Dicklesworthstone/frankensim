use super::*;
use super::super::{drum_with_material,Stroke,acoustics,snare};
use fs_exec::CancelGate;

fn text()->String {format!("{HEADER}\ninitial,relaxed\nband_hz,20,2000\nhead,batter,4000000000\nbranch,batter,1000000000,0.001\nhead,resonant,4000000000\nbranch,resonant,600000000,0.002\n")}
fn specimen()->drum_spec::Spec {
    let mut spec=drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()};
    for head in &mut spec.heads {head.damping_ratio=0.0;}
    spec
}
fn energy(system:&Mechanics)->f64 {
    match system {
        Mechanics::Reference(s)=>s.stored_energy_j(),Mechanics::Nonlinear(s)=>s.stored_energy_j(),
        Mechanics::Substepped(s)=>s.stored_energy_j(),Mechanics::Driven {inner,..}=>energy(inner),
        Mechanics::Prepared(_)=>panic!("expected hereditary-capable mechanics"),
    }
}

#[test]
fn complete_material_inputs_refuse_duplicates_missing_data_and_double_counted_loss() {
    let valid=text();let material=Spec::read(&valid).unwrap();material.admit(&specimen()).unwrap();
    assert!(material.admit(&drum_spec::Spec::reference()).is_err());
    for bad in [valid.replace("head,batter,4000000000\n",""),valid.replace("initial,relaxed\n",""),
        valid.replace("band_hz,20,2000","band_hz,2000,20"),valid.replace("0.001","NaN"),
        valid.replace("1000000000","-1"),format!("{valid}head,resonant,4000000000\n"),
        format!("{valid}initial,unrelaxed\n"),format!("{valid}invented,1\n"),
        format!("{valid}{}","branch,batter,1,0.1\n".repeat(8))] {assert!(Spec::read(&bad).is_err());}
    let mut changed=specimen();changed.heads[0].young_pa*=2.0;assert!(material.admit(&changed).is_err());
    let mut args=vec!["drum".into(),"--head-relaxation".into(),"film.fshr".into(),"128".into()];
    assert_eq!(option(&mut args).unwrap().as_deref(),Some("film.fshr"));assert_eq!(args,["drum","128"]);
    assert!(option(&mut vec!["--head-relaxation".into(),"--analytic-newton".into()]).is_err());
    for name in ["drum","drum-mic","drum-stretch-wav","snare","snare-off-mic"] {admit_command(true,name).unwrap();}
    for name in ["splash","drum-modal","drum-modal-mic"] {assert!(admit_command(true,name).is_err());}
}

#[test]
fn bending_energy_factor_excludes_installed_tension_and_matches_spatial_work() {
    let spec=specimen();let (films,modes)=spec.prepare(2e-6,false).unwrap();let film=&films[0];
    let (matrix,l)=bending_factor(film,&modes[0]).unwrap();let n=modes[0].len();
    let q:Vec<_>=(0..n).map(|i|0.001/(i+1) as f64).collect();
    let reduced=0.5*(0..n).map(|i|q[i]*(0..n).map(|j|matrix[i*n+j]*q[j]).sum::<f64>()).sum::<f64>();
    let factored=0.5*(0..n).map(|i|(0..n).map(|j|l.l(j,i)*q[j]).sum::<f64>().powi(2)).sum::<f64>();
    let rim:Vec<_>=(0..film.mesh.nodes.len()).filter(|&i|film.model.dof_map[3*i].is_none()).collect();
    let model=fs_plate::assemble(&film.mesh,&film.section,&rim,&[],&fs_plate::AssemblyOptions {
        pretension:0.0,support:fs_plate::EdgeSupport::SimplySupported}).unwrap();
    let mut u=vec![0.0;model.free];for (mode,q) in modes[0].iter().zip(&q) {for (u,phi) in u.iter_mut().zip(&mode.phi) {*u+=q*phi;}}
    let mut force=vec![0.0;model.free];model.k.spmv(&u,&mut force);
    let spatial=0.5*u.iter().zip(force).map(|(u,f)|u*f).sum::<f64>();
    assert!(spatial>0.0 && (reduced-spatial).abs()<1e-10*spatial && (factored-spatial).abs()<1e-10*spatial);
    let mut different=spec;different.heads[0].tension_n_m*=10.0;
    let higher=different.head(0).unwrap();let (same,_)=bending_factor(&higher,&modes[0]).unwrap();
    assert_eq!(matrix,same,"material bending must not absorb the prestress operator");
    let total=0.5*modes[0].iter().zip(q).map(|(mode,q)|mode.lambda*q*q).sum::<f64>();assert!(total>spatial);
    // The complex constitutive response is the fs-material law, at every
    // coupling entry, not a constant loss ratio chosen at one resonance.
    let material=Spec::read(&text()).unwrap();let law=&material.heads[0];
    for omega in [30.0,300.0,3000.0] {
        let (storage,loss)=law.modulus(omega);let t=omega*law.terms[0].1;
        assert!((loss/law.e_inf-law.terms[0].0/law.e_inf*t/(1.0+t*t)).abs()<1e-14);
        assert!(storage>law.e_inf && loss>0.0);
    }
}

#[test]
fn real_stick_contact_drives_material_memory_and_changes_head_motion_not_observer_gain() {
    let material=Spec::read(&text()).unwrap();
    let build=|law|drum_with_material(256,2e-6,false,false,None,true,
        Stroke {speed_m_s:2.0,position_m:Some([0.06,0.01])},false,None,Some(specimen()),None,&[],0.0,None,false,law).unwrap();
    let mut elastic=build(None);let mut relaxing=build(Some(&material));let prefix=elastic.system.state().len();
    assert_eq!(elastic.system.state(),&relaxing.system.state()[..prefix]);
    assert_eq!(elastic.force,relaxing.force);assert_eq!(elastic.observer_a,relaxing.observer_a);
    assert_eq!(elastic.pressure.as_ref().unwrap().areas,relaxing.pressure.as_ref().unwrap().areas);
    let initial=energy(&relaxing.system);let mut loss=0.0;let mut changed=0.0_f64;let mut memory_peak=0.0_f64;
    elastic.system=elastic.system.into_analytic_nonlinear().unwrap();relaxing.system=relaxing.system.into_analytic_nonlinear().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..256 {
        let frame=relaxing.system.step(&relaxing.force,&gate).unwrap();elastic.system.step(&elastic.force,&gate).unwrap();
        loss+=frame.dissipated_energy_j;
        for (a,b) in elastic.system.state().iter().zip(relaxing.system.state()) {changed=changed.max((a-b).abs());}
        let memory=observation(&relaxing.system);memory_peak=memory_peak.max(memory.stored_energy_j);
        assert!(memory.dissipated_power_w>=0.0);assert!(frame.balance_residual_j.abs()<1e-7);
    }
    assert!(changed>0.0 && memory_peak>0.0 && loss>0.0);
    assert!((energy(&relaxing.system)-initial+loss).abs()<1e-6);
}

#[test]
fn complete_snare_and_air_keep_their_addresses_before_material_history() {
    let material=Spec::read(&text()).unwrap();
    let build=|law|drum_with_material(2,acoustics::MECHANICAL_DT,true,false,Some(snare::SnareSet::reference(false)),true,
        Stroke::default(),true,None,Some(specimen()),Some(Stroke {speed_m_s:0.0,position_m:Some([-0.05,0.02])}),
        &[],25.0,None,false,law).unwrap();
    let original=build(None);let hereditary=build(Some(&material));
    assert_eq!(original.force,hereditary.force);assert_eq!(original.observer_a,hereditary.observer_a);
    assert_eq!(original.system.state(),&hereditary.system.state()[..original.system.state().len()]);
    assert_eq!(original.acoustics.as_ref().unwrap().state_modes(),hereditary.acoustics.as_ref().unwrap().state_modes());
    assert!(hereditary.force.len()>160);
    assert_eq!(hereditary.system.state().len()-original.system.state().len(),observation(&hereditary.system).states);
    let mut excessive=material.clone();for law in &mut excessive.heads {law.terms=vec![(1e16,0.001)];}
    assert!(drum_with_material(1,2e-6,false,false,None,false,Stroke::default(),false,None,Some(specimen()),
        None,&[],0.0,None,false,Some(&excessive)).is_err());
}

#[test]
fn explicit_elastic_material_and_invalid_numerical_image_do_not_hide_changes() {
    let mut material=Spec::read(&text()).unwrap();for law in &mut material.heads {law.terms.clear();}
    let build=|law|drum_with_material(8,2e-6,false,false,None,false,Stroke::default(),false,None,
        Some(specimen()),None,&[],0.0,None,false,law).unwrap();
    let mut original=build(None);let mut selected=build(Some(&material));let gate=CancelGate::new_clock_free();
    for _ in 0..8 {original.system.step(&original.force,&gate).unwrap();selected.system.step(&selected.force,&gate).unwrap();
        assert_eq!(original.system.state(),selected.system.state());}
    assert_eq!(observation(&selected.system).states,0);
    assert!(drum_with_material(1,2e-6,false,true,None,false,Stroke::default(),false,None,Some(specimen()),
        None,&[],0.0,None,false,Some(&material)).is_err());
}
