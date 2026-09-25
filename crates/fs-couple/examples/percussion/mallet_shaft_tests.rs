use super::*;
use crate::{shaft_playing, Stroke};
use fs_couple::render::plate::impact::{ImpactBody, ImpactSystem};
use fs_couple::render::plate::impact::striker::RadiusStation;
use fs_exec::CancelGate;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;

fn shaft() -> FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.,maximum_modes:17},0.001).unwrap()
}
fn spec() -> Spec {Spec::parse(include_str!("estimated-loaded-mallet.fsmallet")).unwrap()}
fn selection() -> shaft_playing::Selection {shaft_playing::Selection{first:Some(shaft()),second:None}}
fn stroke() -> Stroke {Stroke{speed_m_s:0.5,position_m:Some([0.06,0.01])}}
fn prepared_tip() -> (FeltStriker, Built) {
    let spec=spec();let tips=crate::mallets::Selection{first:Some(spec.clone()),second:None};
    let mut built=selection().build_with_mallets(2,2,stroke(),None,2e-6,&tips).unwrap();
    let tip=spec.from_rows(&vec![vec![1.0/0.1_f64.sqrt()];4],0,built.total,0.5).unwrap();
    (spec.bind_shaft(tip,0,&mut built).unwrap(),built)
}

#[test]
fn physical_head_input_never_reinterprets_a_v1_effective_mass() {
    let v1=Spec::parse(include_str!("estimated-felt-mallet.fsmallet")).unwrap();
    assert!(v1.admit_shaft(false).is_ok());assert!(v1.admit_shaft(true).is_err());
    let v2=spec();assert!(v2.admit_shaft(true).is_ok());assert!(v2.admit_shaft(false).is_err());
    let text=include_str!("estimated-loaded-mallet.fsmallet");
    for bad in [text.replace("attachment,0.0000012,0",""),
        text.replace("attachment,0.0000012,0","attachment,-0.1,0"),
        text.replace("attachment,0.0000012,0","attachment,0.0000012,NaN"),
        text.replace("attachment,0.0000012,0","attachment,0.0000012,4"),
        text.replace("-v2","-v1"),format!("{text}\nattachment,0.0000012,0")] {
        assert!(Spec::parse(&bad).is_err());
    }
    assert!(v2.attach_inertia(&v2.attach_inertia(&shaft()).unwrap()).is_err());
}

#[test]
fn finite_face_preserves_initial_gaps_and_applies_force_plus_moment_transposes() {
    let spec=spec();let (tip,built)=prepared_tip();let p=built.ports[0].as_ref().unwrap();
    let n=built.total;let b=p.tip_row(n).unwrap();let r=p.tip_slope_row(n).unwrap();
    let mut state=vec![0.;2*n];state[0]=tip.body.initial[0].displacement_m_sqrt_kg;
    for i in 0..n {state[2*i+1]=0.003*(i+1) as f64;}
    let o=p.observe(&state).unwrap();let a=spec.radius_m/std::f64::consts::SQRT_2;
    let loads=[1.0,0.4,0.2,0.7];let levers=[a,0.,-a,0.];let mut power=0.;
    for (i,pad) in tip.pads.iter().enumerate() {
        assert!((pad.precompression_m+pad.weights[0]*state[0]+spec.jaw.initial_gap_m).abs()<1e-18);
        assert_eq!(pad.weights[1],-1.0/0.1_f64.sqrt());
        for j in [0].into_iter().chain(p.elastic_start()..p.elastic_start()+p.elastic_modes()) {
            assert!((pad.weights[j]-b[j]-levers[i]*r[j]).abs()<1e-12);
        }
        assert_eq!(pad.area_m2,std::f64::consts::PI*spec.radius_m.powi(2)/4.);
        assert_eq!(pad.creep[0].stiffness_n_m,750.);
        power+=loads[i]*(0..n).map(|j|pad.weights[j]*state[2*j+1]).sum::<f64>();
    }
    let target_speed=state[3]/0.1_f64.sqrt();
    let force=loads.iter().sum::<f64>();let moment=loads.iter().zip(levers).map(|(f,x)|f*x).sum::<f64>();
    assert!((power-force*(o.tip_velocity_m_s-target_speed)-moment*o.tip_angular_velocity_rad_s).abs()<1e-12);
    assert!(moment>0.);assert!(built.bodies[0].is_none(),"no duplicate rigid head body remains");
    assert_eq!(built.elastic.len(),1);
}

#[test]
fn loaded_felt_contact_excites_flexure_and_preserves_histories_on_rejection() {
    let make=|| {
        let (tip,mut built)=prepared_tip();let ports=built.ports[0].take().unwrap();let n=built.total;
        let target=ImpactBody::free_mass(0.1,0.0,0.0).unwrap().0;
        let system=ImpactSystem::new(vec![tip.body,target,built.elastic.remove(0)],vec![],tip.pads,
            vec![],crate::config(512,2e-6)).unwrap().prepare_analytic().unwrap();
        (system,ports,n)
    };
    let (mut a,p,n)=make();let (mut clean,_,_)=make();let force=vec![0.;n];
    let gate=CancelGate::new_clock_free();let initial=a.stored_energy_j();let mut loss=0.;
    let mut flexure=0.0_f64;
    for tick in 0..512 {
        if tick==180 {
            let before=a.state().to_vec();let h=a.felt_history(0).unwrap();
            let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(a.step(&force,&cancel).is_err());
            let mut bad=force.clone();bad[0]=1e9;assert!(a.step(&bad,&gate).is_err());
            assert_eq!(a.state(),before);assert_eq!(a.felt_history(0).unwrap(),h);
        }
        let f=a.step(&force,&gate).unwrap();clean.step(&force,&gate).unwrap();
        assert_eq!(a.state(),clean.state());loss+=f.dissipated_energy_j;
        assert_eq!(f.supplied_work_j,0.);assert!(f.balance_residual_j.abs()<1e-7);
        flexure=flexure.max(p.observe(a.state()).unwrap().flexural_energy_j);
    }
    assert!(flexure>1e-12 && loss>0.);
    assert!((a.stored_energy_j()+loss-initial).abs()<1e-6);
    assert!(a.felt_history(0).unwrap().eps_max>0.);
}

#[test]
fn both_loaded_heads_compose_with_real_drum_snare_cavity_and_shell_sources() {
    let tips=crate::mallets::Selection{first:Some(spec()),second:Some(spec())};
    let shafts=shaft_playing::Selection{first:Some(shaft()),second:Some(shaft())};
    let second=Some(Stroke{speed_m_s:0.3,position_m:Some([-0.05,0.02])});
    let drum=crate::drum_spec::Spec{radial_intervals:2,azimuths:8,..crate::drum_spec::Spec::reference()};
    let wire=crate::snare::SnareSet{strands:2,modes_per_strand:2,contact_cells:4,
        ..crate::snare::SnareSet::reference(false)};
    let d=crate::drum_with_shafts(32,2e-6,true,false,Some(wire),true,stroke(),true,None,
        Some(drum),second,&[],20.,None,false,None,&tips,&shafts).unwrap();
    let mut shell=crate::specimen::Specimen::reference();shell.azimuths=8;
    let c=crate::splash_with_shafts(32,2e-6,true,stroke(),Some(shell),&[],second,None,&tips,&shafts).unwrap();
    for e in [d,c] {
        let sources=e.acoustics.as_ref().unwrap().state_modes();
        let last_rigid=e.second_stick.unwrap().coordinate;
        assert!(sources.iter().all(|i|*i>0 && *i<last_rigid));
        for p in e.flexible_sticks.iter().flatten() {
            let o=p.observe(e.system.state()).unwrap();
            assert!((o.tip_displacement_m+spec().jaw.initial_gap_m).abs()<1e-17);
            assert_eq!(o.flexural_energy_j,0.);
            for i in p.elastic_start()..p.elastic_start()+p.elastic_modes() {
                assert!(!sources.contains(&i));
                if let Some(v)=&e.pressure {assert_eq!(v.areas[i],0.);}
            }
        }
        let crate::Mechanics::Reference(system)=&e.system else{panic!()};
        for i in 0..8 {assert!(system.felt_history(i).is_some());}
    }
}
