use super::*;
use super::super::super::{board, geometry};

fn courses() -> Vec<Course> {
    let all = geometry::demonstration_scale().unwrap();
    vec![Course { strike_fraction: 0.5, ..all[48] }]
}
fn spec(courses: &[Course], span: bool) -> Specification {
    Specification { footprints: courses.iter().map(|c| (c.midi,
        if span { Footprint::Span { length_m: 0.08*c.length_m, sites: 4 } }
        else { Footprint::Point })).collect() }
}
fn bank(span: bool, damped: bool) -> Bank {
    let c=courses();
    Bank::new_with_hammer_footprints(&c,&board::demonstration(),192_000,21_600.0,12,damped,&spec(&c,span)).unwrap()
}

#[test]
fn footprint_input_is_complete_physical_and_quadrature_conserves_area() {
    let c=courses();
    let valid=format!("{HEADER}\nspan,69,0.02,4\n");
    Specification::read(&valid,&c).unwrap();
    for text in [String::new(),format!("{HEADER}\n"),valid.replace("0.02","NaN"),
        valid.replace("0.02","100"),valid.replace("0.02","-1"),valid.replace(",4",",3"),
        valid.replace(",69,",",68,"),format!("{valid}point,69\n"),
        valid.replace("span,69,0.02,4","point,69,0.02,4")] {
        assert!(Specification::read(&text,&c).is_err(),"accepted {text}");
    }
    for sites in [2,4] {
        let q=quadrature(sites);
        assert!(q.iter().all(|(x,w)| x.abs()<1.0 && *w>0.0));
        for power in 0..2*sites {
            let actual=q.iter().map(|(x,w)|w*x.powi(power as i32)).sum::<f64>();
            let expected=if power%2==0 {1.0/(power+1) as f64}else{0.0};
            assert!((actual-expected).abs()<2e-15);
        }
    }
}

#[test]
fn authored_profile_requires_complete_ordered_geometry_and_conserved_area() {
    let c=courses();
    let one=format!("{HEADER}\nprofile,69\nsite,69,0,0,0.008,1\n");
    Specification::read(&one,&c).unwrap();
    let four=format!("{HEADER}\nprofile,69\nsite,69,-0.003,0.0003,0.006,0.15\nsite,69,-0.001,0,0.008,0.35\nsite,69,0.001,0,0.008,0.35\nsite,69,0.003,0.0003,0.006,0.15\n");
    Specification::read(&four,&c).unwrap();
    for invalid in [
        one.replace("profile,69\n",""),format!("{HEADER}\nprofile,69\n"),
        format!("{one}point,69\n"),format!("{four}site,69,0.004,0,0.008,0.1\n"),
        one.replace("profile,69","point,69"),one.replace("site,69","site,68"),
        one.replace("site,69,0,0,0.008,1","site,69,100,0,0.008,1"),
        one.replace("site,69,0,0,0.008,1","site,69,0,-0.001,0.008,1"),
        one.replace("0.008","0"),one.replace("0.008","NaN"),
        one.replace(",0,0.008",",inf,0.008"),
        one.replace("0.008,1","0.008,0"),one.replace("0.008,1","0.008,0.9"),
        four.replace("site,69,0.001","site,69,-0.001"),
        four.replace("site,69,0.001","site,69,-0.002"),
    ] {assert!(Specification::read(&invalid,&c).is_err(),"accepted {invalid}");}
}

#[test]
fn profile_projects_supplied_stations_without_putting_recession_in_string_work() {
    let c=courses();
    let text=format!("{HEADER}\nprofile,69\nsite,69,-0.004,0.0004,0.006,0.3\nsite,69,0.002,0,0.008,0.7\n");
    let selection=Specification::read(&text,&c).unwrap();
    let mut b=Bank::new_with_hammer_footprints(&c,&board::demonstration(),192_000,
        21_600.0,12,true,&selection).unwrap();
    let point=bank(false,true);
    assert_eq!(b.q,point.q);assert_eq!(b.board_basis,point.board_basis);
    assert_eq!(b.physical_board_k,point.physical_board_k);
    assert_eq!(b.contact_strings.len(),2*point.contact_strings.len());
    for (i,(&offset,(&recession,&thickness))) in [-0.004,0.002].iter()
        .zip([0.0004,0.0].iter().zip([0.006,0.008].iter())).enumerate() {
        assert_eq!(b.contact_recession_m(i),recession);
        assert_eq!(b.contact_felt_thickness_m(i),Some(thickness));
        // Static face geometry cannot do work on a stationary string.
        assert_eq!(b.contact_position(i,&b.q),0.0);
        let si=b.contact_strings[i];let k=b.strings[si].modes.start+1;
        b.q[k]=1e-5;
        let station=c[0].strike_fraction+offset/c[0].length_m;
        let expected=det::sin(2.0*std::f64::consts::PI*station)
            /det::sqrt(c[0].modal_mass_kg())*b.q[k];
        assert!((b.contact_position(i,&b.q)-expected).abs()<1e-15);
        b.q[k]=0.0;
    }
    // Arbitrary signed forces check the actual moving-bridge work row, beyond
    // a getter or a duplicate projection formula.
    let force:Vec<_>=(0..b.contact_strings.len()).map(|i|0.2*(i as f64+0.5).sin()).collect();
    for _ in 0..12 {
        let before=b.energy();b.predict();b.finish(&force);
        let work:f64=force.iter().enumerate().map(|(i,f)|
            f*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum();
        assert!((b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work).abs()<1e-11);
        b.commit();
    }
}

#[test]
fn finite_sites_keep_every_mode_and_do_not_average_away_a_nodal_contact() {
    let mut spatial=bank(true,false);let point=bank(false,false);
    assert_eq!(spatial.q.len(),point.q.len());assert_eq!(spatial.board_basis,point.board_basis);
    assert_eq!(spatial.physical_board_k,point.physical_board_k);
    assert_eq!(spatial.strings.len(),point.strings.len());
    assert_eq!(spatial.contact_strings.len(),4*point.contact_strings.len());
    let si=spatial.contact_strings[0];let s=&spatial.strings[si];
    spatial.q[s.modes.start+1]=1e-5; // second sine mode: node at face centre
    let a=spatial.contact_position(0,&spatial.q);let b=spatial.contact_position(3,&spatial.q);
    assert!(a*b<0.0 && a.abs()>1e-7);
    assert!(point.contact_position(0,&spatial.q).abs()<1e-16);
    let fraction:f64=(0..spatial.contact_strings.len()).filter(|i|spatial.contact_strings[*i]==si)
        .map(|i|spatial.contact_area_fraction(i)).sum();
    assert!((fraction-1.0).abs()<1e-15);
    assert!(spatial.strings.iter().filter(|s|s.contact.is_none()).all(|s|
        !spatial.contact_strings.iter().any(|i|spatial.strings[*i].modes==s.modes)));
}

#[test]
fn every_site_compliance_matches_the_actual_coupled_force_response() {
    let mut b=bank(true,true);let nc=b.contact_strings.len();let mut force=vec![0.0;nc];
    for c in 0..nc {
        force[c]=1.0;b.predict();b.finish(&force);
        for i in 0..nc {
            let measured=b.contact_position(i,&b.next_q)-b.free_contact[i];
            assert!((measured-b.contact_compliance[i*nc+c]).abs()<1e-14);
            assert_eq!(b.contact_compliance[i*nc+c],b.contact_compliance[c*nc+i]);
        }
        force[c]=0.0;
    }
    assert!(b.contact_compliance[1].abs()>1e-10,"same-string cross-site compliance must not be omitted");
}

#[test]
fn arbitrary_site_forces_close_the_original_string_board_work_balance() {
    for damped in [false,true] {
        let mut b=bank(true,damped);let nc=b.contact_strings.len();
        for k in 0..b.q.len() {b.q[k]=1e-7*(k as f64).sin();b.v[k]=0.001*(k as f64).cos();}
        for tick in 0..100 {
            let f:Vec<_>=(0..nc).map(|i|0.1*((i+tick) as f64).sin()).collect();
            let before=b.energy();b.predict();b.finish(&f);
            let work:f64=(0..nc).map(|c|f[c]*(b.contact_position(c,&b.next_q)-b.contact_position(c,&b.q))).sum();
            let defect=b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work;
            assert!(defect.abs()<1e-11,"{defect:e}");b.commit();
        }
    }
}

#[test]
fn explicit_point_selection_preserves_legacy_motion_and_failed_preparation_is_atomic() {
    let c=courses();let mut plain=Bank::new(&c,&board::demonstration(),192_000,21_600.0,12,true).unwrap();
    let mut selected=bank(false,true);let force=vec![0.2;plain.contact_strings.len()];
    for _ in 0..20 {
        plain.predict();selected.predict();plain.finish(&force);selected.finish(&force);
        plain.commit();selected.commit();assert_eq!(plain.q,selected.q);assert_eq!(plain.v,selected.v);
    }
    let before=selected.q.clone();let compliance=selected.contact_compliance.clone();
    assert!(selected.configure_hammer_footprints(&spec(&c,true),&c).is_err());
    assert_eq!(selected.q,before);assert_eq!(selected.contact_compliance,compliance);
}
