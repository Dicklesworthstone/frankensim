use super::*;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;
use fs_couple::render::plate::impact::striker::RadiusStation;

fn shaft() -> FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.0,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.0,maximum_modes:17},0.001).unwrap()
}
fn specimen() -> specimen::Specimen {
    let mut s=specimen::Specimen::reference();s.azimuths=8;s
}
fn spec() -> Spec {Spec::parse(include_str!("estimated-hihat.fshh")).unwrap()}
fn strokes() -> [Stroke;2] {
    [Stroke{speed_m_s:0.8,position_m:Some([0.06,0.01])},
     Stroke{speed_m_s:0.6,position_m:Some([-0.06,0.01])}]
}

#[test]
fn flexible_shafts_append_without_shifting_shell_pedal_or_receiver_sources() {
    let shaft=shaft();let [a,b]=strokes();let specimen=specimen();
    let pair=build_with_strikers(&spec(),&specimen,&specimen,a,Some(b),256,2e-6,true,None,
        [Some(&shaft),Some(&shaft)]).unwrap();
    let e=&pair.experiment;let n=e.force.len();let p=pair.flexible_sticks[0].as_ref().unwrap();
    let q=pair.flexible_sticks[1].as_ref().unwrap();let second=e.second_stick.unwrap();
    assert_eq!(p.rigid_coordinate(),0);assert_eq!(q.rigid_coordinate(),second.coordinate);
    assert_eq!(second.coordinate,pair.pedal.coordinate+1);
    assert_eq!(p.elastic_start(),second.coordinate+1);
    assert_eq!(q.elastic_start(),p.elastic_start()+p.elastic_modes());
    assert_eq!(n,q.elastic_start()+q.elastic_modes());
    assert_eq!(e.system.state().len(),2*n+12);
    let sources=e.acoustics.as_ref().unwrap().state_modes();
    assert_eq!(sources,pair.upper_modes.clone().chain(pair.lower_modes.clone()).collect::<Vec<_>>());
    for port in [p,q] {
        assert!(!sources.contains(&port.rigid_coordinate()));
        for k in port.elastic_start()..port.elastic_start()+port.elastic_modes(){assert!(!sources.contains(&k));}
        let o=port.observe(e.system.state()).unwrap();
        assert_eq!(o.flexural_energy_j,0.0);assert!((o.tip_displacement_m+0.0002).abs()<1e-16);
        assert_ne!(port.tip_row(n).unwrap(),port.hand_row(n).unwrap());
        for row in pair.collision.collocation().chunks(n) {
            assert_eq!(row[port.rigid_coordinate()],0.0);
            assert!(row[port.elastic_start()..port.elastic_start()+port.elastic_modes()].iter().all(|x|*x==0.0));
        }
    }
    assert!((p.observe(e.system.state()).unwrap().tip_velocity_m_s-a.speed_m_s).abs()<1e-14);
    assert!((q.observe(e.system.state()).unwrap().tip_velocity_m_s-b.speed_m_s).abs()<1e-14);
}

#[test]
fn actual_hi_hat_contacts_excite_both_shafts_inside_the_shared_energy_solve() {
    let shaft=shaft();let [a,b]=strokes();let specimen=specimen();
    let pair=build_with_strikers(&spec(),&specimen,&specimen,a,Some(b),256,2e-6,false,None,
        [Some(&shaft),Some(&shaft)]).unwrap();
    let mut e=pair.experiment;
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        fs_couple::render::plate::impact::ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();
    let gate=CancelGate::new_clock_free();let mut bending=[0.0_f64;2];
    let (mut energy,mut loss,mut initial)=(0.0,0.0,None);
    for tick in 0..256 {
        if tick==160 {let before=e.system.state().to_vec();let mut bad=e.force.clone();bad[0]=f64::NAN;
            assert!(e.system.step(&bad,&gate).is_err());assert_eq!(e.system.state(),before);}
        let f=e.system.step(&e.force,&gate).unwrap();
        if initial.is_none(){initial=Some(f.stored_energy_j+f.dissipated_energy_j-f.supplied_work_j);}
        assert!(f.balance_residual_j.abs()<1e-7);assert_eq!(f.supplied_work_j,0.0);
        loss+=f.dissipated_energy_j;energy=f.stored_energy_j;
        for (i,port) in pair.flexible_sticks.iter().enumerate(){
            bending[i]=bending[i].max(port.as_ref().unwrap().observe(e.system.state()).unwrap().flexural_energy_j);
        }
    }
    assert!(bending.iter().all(|x|*x>1e-12),"both actual contacts must excite retained flexure: {bending:?}");
    assert!((energy+loss-initial.unwrap()).abs()<1e-6);
}

#[test]
fn second_only_flexible_selection_keeps_the_original_first_stick() {
    let shaft=shaft();let [a,b]=strokes();let specimen=specimen();
    let pair=build_with_strikers(&spec(),&specimen,&specimen,a,Some(b),1,2e-6,false,None,
        [None,Some(&shaft)]).unwrap();
    assert!(pair.flexible_sticks[0].is_none());
    let p=pair.flexible_sticks[1].as_ref().unwrap();
    assert_eq!(p.elastic_start(),p.rigid_coordinate()+1);
    assert!((pair.experiment.system.state()[1]*pair.experiment.stick_weight-a.speed_m_s).abs()<1e-14);
    assert!(build_with_strikers(&spec(),&specimen,&specimen,a,None,1,2e-6,false,None,
        [None,Some(&shaft)]).is_err());
}
