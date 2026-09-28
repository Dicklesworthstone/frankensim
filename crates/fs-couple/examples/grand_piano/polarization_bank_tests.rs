//! Two physical transverse coordinates per span, not duplicated struck voices.
use super::*;
fn course()->Course {
    let c=super::super::geometry::demonstration_scale().unwrap()[48];
    Course {unison:1,duplex_length_m:0.,detune_cents:0.,..c}
}
fn board()->Vec<BoardMode> {
    [170.,310.].into_iter().enumerate().map(|(i,f)| {
        let mut bridge=[0.;88];bridge[48]=if i==0 {0.08}else{-0.04};
        BoardMode {frequency_hz:f,damping_ratio:0.01,bridge,volume:0.1}
    }).collect()
}
fn pair(side:f64,damping:bool)->Bank {
    Bank::new_with_transverse_bridge(&[course()],&board(),192_000,21_600.,8,damping,
        Some(&[vec![side,0.5*side]])).unwrap()
}
fn step(b:&mut Bank,force:f64) {
    let before=b.energy();b.begin_string_stretching_step();
    let mut converged=false;
    for _ in 0..string_stretching::MAX_TRIALS {
        b.predict();b.finish(&[force]);
        if b.correct_string_stretching_step().unwrap() {converged=true;break;}
    }
    assert!(converged);
    let work=force*(b.contact_position(0,&b.next_q)-b.contact_position(0,&b.q));
    let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
    assert!(defect.abs()<1e-10,"two-polarization work defect: {defect:e}");
    b.commit();
}
#[test]
fn second_direction_is_reciprocal_and_not_a_second_hammer_or_retuned_string() {
    let mut b=pair(0.06,true);
    assert!(b.has_secondary_polarization());assert_eq!(b.contact_strings.len(),1);
    assert_eq!(b.strings.len(),2);assert_eq!(b.groups.len(),2);
    let a=&b.strings[0];let c=&b.strings[1];
    assert_eq!(a.member,c.member);assert_eq!(a.course,c.course);
    assert_eq!(a.duplex,c.duplex);assert_eq!((a.polarization,c.polarization),(0,1));
    assert!(c.contact.is_none());assert!(!c.duplex);
    for (i,j) in a.modes.clone().zip(c.modes.clone()) {
        assert_eq!(b.modes[i].omega,b.modes[j].omega);
        assert_eq!(b.modes[i].beta,b.modes[j].beta);
    }
    let second=c.modes.clone();
    for n in 0..1000 {step(&mut b,if n<100 {0.2}else{0.});}
    assert!(b.q[second].iter().any(|v|v.abs()>1e-12),"unstruck transverse motion must receive bridge work");
}
#[test]
fn mass_completion_retains_the_full_vector_kinetic_energy_and_port_work() {
    let b=pair(0.06,false);let n=b.modes.len();let r=b.board_count;let c=course();
    let v:Vec<_>=(0..n+r).map(|i|0.002*(i as f64+0.3).sin()).collect();
    let bare:Vec<_>=(0..r).map(|i|(0..r).map(|j|b.board_basis[i*r+j]*v[n+j]).sum::<f64>()).collect();
    let mut kinetic=0.5*bare.iter().map(|v|v*v).sum::<f64>();
    for s in &b.strings {
        let speed=s.bridge.iter().zip(&v[n..]).map(|(g,v)|g*v).sum::<f64>();
        let mut cross=0.;
        for k in s.modes.clone() {
            let relative=v[k]-b.modes[k].beta*speed;
            kinetic+=0.5*relative*relative;cross+=b.modes[k].beta*relative;
        }
        kinetic+=cross*speed+0.5*(c.linear_density_kg_m*c.length_m/3.)*speed*speed;
    }
    let completed=0.5*v.iter().map(|v|v*v).sum::<f64>();
    assert!((kinetic-completed).abs()<1e-12*completed);
    // Check the actual scalar contact compliance after both endpoint mass forms.
    let mut b=b;b.predict();b.finish(&[1.]);
    assert!((b.contact_position(0,&b.next_q)-b.free_contact[0]-b.contact_compliance[0]).abs()<1e-14);
}
#[test]
fn a_zero_transverse_port_keeps_the_original_primary_trajectory_and_observer() {
    let mut a=Bank::new(&[course()],&board(),192_000,21_600.,8,true).unwrap();
    let mut b=pair(0.,true);let na=a.modes.len();let nb=b.modes.len();
    assert_eq!(a.contact_compliance,b.contact_compliance);
    for frame in 0..300 {
        let f=if frame<50 {0.2}else{0.};step(&mut a,f);step(&mut b,f);
        assert_eq!(&a.q[..na],&b.q[..na]);assert_eq!(&a.v[..na],&b.v[..na]);
        assert_eq!(&a.q[na..],&b.q[nb..]);assert_eq!(&a.v[na..],&b.v[nb..]);
        assert_eq!(a.volume_velocity(),b.volume_velocity());
        assert!(b.q[na..nb].iter().chain(&b.v[na..nb]).all(|v|*v==0.));
    }
}
#[test]
fn orthogonal_motion_shares_one_extension_tension_and_cross_polarization_energy() {
    let mut b=pair(0.06,false);let c=course();
    let spec=string_stretching::Specification::read(
        "frankensim-piano-string-stretching-v1\nstretch,69,100000,0.2\n",&[c]).unwrap();
    b.configure_string_stretching(&[c],&spec).unwrap();
    let a=b.strings[0].modes.start;let d=b.strings[1].modes.start;
    b.q[a]=1e-5;
    let e=b.string_stretching_observation(0).unwrap().stretching_energy_j;
    b.q[d]=1e-5;
    let together=b.string_stretching_observation(0).unwrap();
    let other=b.string_stretching_observation(1).unwrap();
    assert!((together.stretching_energy_j-4.*e).abs()<1e-13*e);
    assert_eq!(together.tension_n,other.tension_n);
    assert_eq!(together.stretching_energy_j,other.stretching_energy_j);
    // Both moving bridge terms participate in the exact discrete-gradient work.
    for k in 0..b.q.len() {b.q[k]=1e-7*(0.7*k as f64).sin();b.v[k]=1e-4*(0.3*k as f64).cos();}
    for _ in 0..100 {step(&mut b,0.1);}
}
#[test]
fn complete_unison_duplex_identity_is_bounded_without_missing_polarizations() {
    let c=Course {unison:3,duplex_length_m:0.12,detune_cents:1.,..course()};
    let b=Bank::new_with_transverse_bridge(&[c],&board(),192_000,21_600.,8,false,
        Some(&[vec![0.03,-0.01]])).unwrap();
    assert_eq!(b.strings.len(),12);assert_eq!(b.contact_strings.len(),3);
    for member in 0..3 {for duplex in [false,true] {
        let pair:Vec<_>=b.strings.iter().filter(|s|s.member==member&&s.duplex==duplex).collect();
        assert_eq!(pair.len(),2);assert_eq!(pair[0].modes.len(),pair[1].modes.len());
    }}
    for wrong in [vec![],vec![vec![0.]],vec![vec![f64::NAN,0.]],vec![vec![0.;2];2]] {
        assert!(Bank::new_with_transverse_bridge(&[c],&board(),192_000,21_600.,8,false,Some(&wrong)).is_err());
    }
    assert!(!Bank::new(&[c],&board(),192_000,21_600.,8,false).unwrap().has_secondary_polarization());
}
