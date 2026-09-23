use super::*;
use super::super::{radiation,relaxation,hammer_footprint};
use super::super::super::{board,geometry};
use fs_material::{WoolFelt,Uniaxial};
use fs_material::visco::GeneralizedMaxwell;

fn piano(spans:bool)->Instrument {
    let scale=geometry::demonstration_scale().unwrap();let courses=vec![scale[48],scale[51]];
    let mut text=format!("{}\n",hammer_footprint::HEADER);
    for c in &courses {text.push_str(&if spans {format!("span,{},{},4\n",c.midi,0.06*c.length_m)}
        else {format!("point,{}\n",c.midi)});}
    let spec=hammer_footprint::Specification::read(&text,&courses).unwrap();
    Instrument::new_with_footprints(courses,&board::demonstration(),48_000,4,12,true,&spec).unwrap()
}
fn balance(p:&Instrument) {
    assert!((p.accounting.input_work_j-p.accounting.dissipated_j()-p.energy_j()).abs()<1e-7);
}
fn same(a:&Instrument,b:&Instrument) {
    assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    assert_eq!(a.radiation_energy_j(),b.radiation_energy_j());
    assert_eq!(a.accounting.input_work_j,b.accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(),b.accounting.dissipated_j());
    for (a,b) in a.contacts.iter().zip(&b.contacts) {
        assert_eq!(a.overlap,b.overlap);assert_eq!(a.force,b.force);assert_eq!(a.memory,b.memory);
        assert_eq!(a.state.eps_max,b.state.eps_max);assert_eq!(a.state.sig_max,b.state.sig_max);
        assert_eq!(a.enabled,b.enabled);
    }
}
fn force_defect(p:&Instrument)->f64 {
    p.active.iter().map(|&i| {
        let ci=p.bank.strings[p.bank.contact_strings[i]].course;
        let material=&p.creep[i];let old=&p.contacts[i];
        let expected=felt::average(&p.laws[ci],&old.state,old.overlap-material.deformation(&old.memory),
            p.gap[i]-material.free_deformation(&old.memory)-material.compliance()*p.force[i],
            p.courses[ci].felt_thickness_m,p.contact_areas[i]).0;
        (expected-p.force[i]).abs()
    }).fold(0.0_f64,f64::max)
}

#[test]
fn actual_loaded_bank_resolves_stiff_nearby_sites_without_collapsing_the_face() {
    let build=|| {
        let mut c=geometry::demonstration_scale().unwrap()[48];c.felt_thickness_m=1e-4;
        let spec=hammer_footprint::Specification::read(&format!("{}\nspan,69,{},4",
            hammer_footprint::HEADER,1e-6*c.length_m),&[c]).unwrap();
        let mut p=Instrument::new_with_contact_geometry(vec![c],&board::demonstration(),8_000,1,8,false,
            vec![(WoolFelt::new(1e9,0.2,2.,2.5,0.2,0.8).unwrap(),
                GeneralizedMaxwell {e_inf:1.,terms:vec![]})],None,Some(&spec)).unwrap();
        let n=p.contacts.len();assert_eq!(n,12);p.active=(0..n).collect();p.force.fill(0.);
        // Manufacture a local contact equation in the ACTUAL condensed bank.
        // No time step is taken with the manufactured free-displacement RHS.
        for old in &mut p.contacts {old.overlap=0.;}
        let wanted:Vec<_>=(0..n).map(|i|felt::average(&p.laws[0],&p.contacts[i].state,
            0.,1e-6,c.felt_thickness_m,p.contact_areas[i]).0).collect();
        for i in 0..n {p.gap[i]=1e-6+(0..n).map(|j|p.contact_h[i*n+j]*wanted[j]).sum::<f64>();}
        (p,wanted)
    };
    let (mut block,wanted)=build();let (mut scalar,_)=build();scalar.contact_solver=None;
    let q=block.bank.q.clone();let v=block.bank.v.clone();
    block.contact_sweep().unwrap();assert!(force_defect(&block)<1e-5);
    for (f,w) in block.force.iter().zip(wanted) {assert!((f-w).abs()<1e-5);}
    let mut refused=false;for _ in 0..32 {if scalar.contact_sweep().is_err(){refused=true;break;}}
    assert!(refused || force_defect(&scalar)>1e-5,"must expose the old scalar-budget limit");
    assert_eq!(block.bank.q,q);assert_eq!(block.bank.v,v);
    assert!(block.contacts.iter().all(|c|c.state.eps_max==0. && c.memory==relaxation::Memory::default()));
}

#[test]
fn two_hammer_chord_preserves_cross_key_work_and_the_converged_scalar_trajectory() {
    let mut block=piano(true);let mut scalar=piano(true);scalar.contact_solver=None;
    for (key,speed) in [(69,2.0),(72,1.4)] {block.note_on(key,speed).unwrap();scalar.note_on(key,speed).unwrap();}
    let mut trace=vec![0.;block.board_trace_len()];let mut reference=trace.clone();
    let mut peak=0.0_f64;
    for _ in 0..1200 {
        block.step_with_board_trace(&mut trace).unwrap();scalar.step_with_board_trace(&mut reference).unwrap();
        peak=peak.max(block.force.iter().copied().fold(0.0_f64,f64::max));
        for (a,b) in block.bank.q.iter().zip(&scalar.bank.q) {assert!((a-b).abs()<1e-8);}
        for (a,b) in trace.iter().zip(&reference) {assert!((a-b).abs()<1e-5);}
    }
    assert!(peak>0.0);assert!(block.contacts.iter().all(|c|c.state.eps_max>0.));
    assert!(block.accounting.felt_relaxation_loss_j>0.);balance(&block);balance(&scalar);
}

#[test]
fn whole_string_pedal_selection_and_failed_second_block_restore_all_histories() {
    let mut a=piano(true);let mut b=piano(true);
    let load=radiation::Model {ports:a.bank.board_count,poles:vec![radiation::Pole {
        omega:std::f64::consts::TAU*300.,zeta:0.02,coupling:vec![40.;a.bank.board_count],
    }]};
    for p in [&mut a,&mut b] {
        p.configure_radiation(&load).unwrap();p.set_una_corda(true);p.set_sustain(0.3).unwrap();
        p.note_on(69,2.).unwrap();p.note_on(72,1.6).unwrap();
    }
    for _ in 0..400 {a.step().unwrap();b.step().unwrap();}
    assert_eq!(a.contacts.iter().filter(|c|c.enabled).count(),16);
    assert!(a.contacts.iter().any(|c|c.state.eps_max>0.));
    let nc=a.contacts.len();let index=(0..nc).find(|&i|a.bank.strings[a.bank.contact_strings[i]].course==1).unwrap();
    // The first hammer's numerical block runs before the malformed second one.
    let diagonal=a.contact_h[index*nc+index];a.contact_h[index*nc+index]=f64::NAN;
    assert!(a.step().is_err());a.contact_h[index*nc+index]=diagonal;same(&a,&b);
    for _ in 0..100 {assert_eq!(a.step().unwrap(),b.step().unwrap());}
    same(&a,&b);balance(&a);
}

#[test]
fn unselected_and_all_point_hammers_keep_the_original_scalar_path_bit_for_bit() {
    let scale=geometry::demonstration_scale().unwrap();
    let mut original=Instrument::new(vec![scale[48],scale[51]],&board::demonstration(),48_000,4,12,true).unwrap();
    let mut points=piano(false);assert!(original.contact_solver.is_none() && points.contact_solver.is_none());
    for key in [69,72] {original.note_on(key,1.5).unwrap();points.note_on(key,1.5).unwrap();}
    for _ in 0..800 {assert_eq!(original.step().unwrap(),points.step().unwrap());}
    same(&original,&points);
}
