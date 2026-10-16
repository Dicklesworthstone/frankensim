use super::*;
use super::super::{BoardMode,hammer_footprint};
use super::super::super::{board,geometry};
use fs_phs::Storage;

fn course()->Course {Course {unison:1,duplex_length_m:0.,detune_cents:0.,..geometry::demonstration_scale().unwrap()[48]}}
fn specification(courses:&[Course],ea:f64)->Specification {
    let mut text=format!("{HEADER}\n");for c in courses {text.push_str(&format!("stretch,{},{ea},0.2\n",c.midi));}
    Specification::read(&text,courses).unwrap()
}
fn bank(courses:&[Course],ea:f64)->Bank {
    let mut b=Bank::new(courses,&board::demonstration(),192_000,21_600.,8,false).unwrap();
    b.configure_string_stretching(courses,&specification(courses,ea)).unwrap();b
}
fn advance(b:&mut Bank,contacts:&[f64])->Result<usize,&'static str> {
    b.begin_string_stretching_step();
    for i in 0..MAX_TRIALS {
        b.predict();b.finish(contacts);
        if b.correct_string_stretching_step()? {return Ok(i+1);}
    }
    Err("nonlinear string trial budget")
}

#[test]
fn complete_material_admission_never_infers_axial_rigidity_from_tension_or_ei() {
    let c=course();let valid=format!("{HEADER}\n# SI inputs\nstretch,69,100000,0.2");
    Specification::read(&valid,&[c]).unwrap();
    for row in ["", "stretch,69,NaN,0.2", "stretch,69,0,0.2", "stretch,69,-1,0.2",
        "stretch,69,100000,0.31", "stretch,69,100000,0", "stretch,70,100000,0.2",
        "linear,69\nlinear,69", "linear,69,1", "nonsense,69"] {
        assert!(Specification::read(&format!("{HEADER}\n{row}"),&[c]).is_err(),"{row}");
    }
    assert!(Specification::read("linear,69",&[c]).is_err());
    let mut b=bank(&[c],100000.);assert!(b.configure_string_stretching(&[c],&specification(&[c],100000.)).is_err());
    b.q[0]=1e-6;
    assert!(b.stretching.as_ref().unwrap().validate_state(&vec![1.;b.q.len()],&b.strings,&b.modes).is_err());
}

#[test]
fn physical_extension_and_virtual_work_include_the_moving_bridge_chord() {
    let c=course();let mut b=bank(&[c],100000.);let n=b.modes.len();
    for (i,q) in b.q.iter_mut().enumerate(){*q=1e-5*(0.7*i as f64).sin();}
    let p=b.stretching.as_mut().unwrap();let ch=&p.channels[0];let s=&b.strings[0];
    let obs=ch.observe(&b.q,&b.strings,&b.modes);
    let endpoint=s.bridge.iter().zip(&b.q[n..]).map(|(g,q)|g*q).sum::<f64>();
    let mut strain=0.;let cells=4096;
    for i in 0..cells {
        let x=(i as f64+0.5)/cells as f64;
        let mut slope=endpoint/c.length_m;
        for (index,k) in s.modes.clone().enumerate() {
            let relative=b.q[k]-b.modes[k].beta*endpoint;
            slope+=relative*(index+1) as f64*std::f64::consts::PI/c.length_m
                *det::cos((index+1) as f64*std::f64::consts::PI*x)/det::sqrt(c.modal_mass_kg());
        }
        strain+=0.5*slope*slope/cells as f64;
    }
    assert!((strain-obs.additional_strain).abs()<1e-15);
    assert!((obs.stretching_energy_j-0.5*100000.*c.length_m*strain*strain).abs()<1e-15);
    p.evaluate(&b.q,&b.q,&b.strings,&b.modes).unwrap();let force=p.required.clone();
    for k in 0..b.q.len() {
        let mut hi=b.q.clone();let mut lo=b.q.clone();let h=1e-9;hi[k]+=h;lo[k]-=h;
        let derivative=(p.energy(&hi,&b.strings,&b.modes)-p.energy(&lo,&b.strings,&b.modes))/(2.*h);
        assert!((force[k]+derivative).abs()<1e-7*(1.+derivative.abs()));
    }
    let mut end=b.q.clone();for (i,q) in end.iter_mut().enumerate(){*q+=2e-6*(i as f64+0.3).cos();}
    p.evaluate(&b.q,&end,&b.strings,&b.modes).unwrap();
    let work=p.required.iter().zip(b.q.iter().zip(&end)).map(|(f,(a,b))|f*(b-a)).sum::<f64>();
    assert!((p.energy(&end,&b.strings,&b.modes)-obs.stretching_energy_j+work).abs()<1e-14);
    // At zero bridge displacement the existing fs-nlmodal potential is exact.
    let owner=kirchhoff_carrier_string(&KcStringParams{length:c.length_m,tension:c.tension_n,
        lin_density:c.linear_density_kg_m,ea:100000.},s.modes.len()).unwrap();
    let mut q=b.q.clone();q[n..].fill(0.);let mut x=vec![0.;2*n];
    for i in 0..n{x[2*i]=q[i];}
    let linear=owner.omegas.iter().zip(&q).map(|(w,q)|0.5*w*w*q*q).sum::<f64>();
    assert!((owner.hamiltonian(&x)-linear-p.energy(&q,&b.strings,&b.modes)).abs()<1e-13);
}

#[test]
fn generalized_trial_force_uses_the_original_reciprocal_propagator_and_loss_ledger() {
    let c=course();let mut b=bank(&[c],100000.);
    // Isolate the reusable affine propagation, not a replacement physical law.
    let p=b.stretching.as_mut().unwrap();p.channels.clear();
    for (i,f) in p.force.iter_mut().enumerate(){*f=0.01*(i as f64+0.2).cos();}
    let external=p.force.clone();let contact=vec![0.3;b.contact_strings.len()];
    for (i,q) in b.q.iter_mut().enumerate(){*q=1e-6*(i as f64).cos();}b.v.fill(0.02);
    let before=b.energy();b.predict();b.finish(&contact);
    let work=external.iter().zip(b.q.iter().zip(&b.next_q)).map(|(f,(a,b))|f*(b-a)).sum::<f64>()
        +(0..contact.len()).map(|i|contact[i]*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum::<f64>();
    assert!((b.energy_at(&b.next_q,&b.next_v)-before+b.last_modal_loss_j-work).abs()<1e-11);
}

#[test]
fn coupled_nonlinear_string_work_closes_without_replacing_unisons_duplexes_or_contacts() {
    let c=geometry::demonstration_scale().unwrap()[48];let mut b=bank(&[c],150000.);
    let shapes=b.modes.iter().map(|m|(m.omega,m.beta,m.hammer_shape)).collect::<Vec<_>>();
    let initial=b.energy();let mut work=0.;let mut loss=0.;
    let rest=b.string_stretching_observation(0).unwrap().tension_n;let mut max_tension=rest;
    for tick in 0..600 {
        let force=vec![if tick<200 {0.5}else{0.};b.contact_strings.len()];
        advance(&mut b,&force).unwrap();
        work+=(0..force.len()).map(|i|force[i]*(b.contact_position(i,&b.next_q)-b.contact_position(i,&b.q))).sum::<f64>();
        loss+=b.last_modal_loss_j;b.commit();
        max_tension=max_tension.max(b.string_stretching_observation(0).unwrap().tension_n);
        assert!((b.energy()-initial+loss-work).abs()<1e-8);
    }
    assert!(max_tension>rest && b.strings.len()==2*c.unison);
    assert_eq!(shapes,b.modes.iter().map(|m|(m.omega,m.beta,m.hammer_shape)).collect::<Vec<_>>());
}

#[test]
fn finite_amplitude_period_converges_to_continuous_duffing_not_a_retuned_frequency() {
    fn period(rate:u32,amplitude:f64)->(f64,f64) {
        let c=course();let board=[BoardMode {frequency_hz:100.,damping_ratio:0.,bridge:[0.;88],volume:0.}];
        let mut b=Bank::new(&[c],&board,rate,3000.,1,false).unwrap();
        b.configure_string_stretching(&[c],&specification(&[c],300000.)).unwrap();
        let ch=&b.stretching.as_ref().unwrap().channels[0];let beta=ch.coefficient*ch.diagonal[0].powi(2);
        let a=amplitude*det::sqrt(c.modal_mass_kg());b.q[0]=a;
        let w=b.modes[0].omega;let mut integral=0.;let cells=4096;
        for i in 0..cells {let theta=(i as f64+0.5)*std::f64::consts::FRAC_PI_2/cells as f64;
            integral+=1./det::sqrt(w*w+0.5*beta*a*a*(1.+det::sin(theta).powi(2)));}
        let exact=4.*std::f64::consts::FRAC_PI_2/cells as f64*integral;
        let initial=b.energy();let mut previous=a;
        for tick in 1..10000 {
            advance(&mut b,&[0.]).unwrap();b.commit();
            assert!((b.energy()-initial).abs()<1e-9);
            if previous>0. && b.q[0]<=0. {
                let crossing=(tick as f64-1.+previous/(previous-b.q[0]))/f64::from(rate);
                return (4.*crossing,exact);
            }
            previous=b.q[0];
        }
        panic!("no zero crossing");
    }
    let weak=period(192000,1e-5);let coarse=period(48000,0.004);let fine=period(96000,0.004);
    assert!(fine.0<weak.0);
    assert!((fine.0-fine.1).abs()<0.5*(coarse.0-coarse.1).abs());
    assert!((fine.0/fine.1-1.).abs()<1e-3);
}

#[test]
fn all_linear_selection_preserves_the_point_and_finite_face_images_bit_for_bit() {
    for spatial in [false,true] {
        let c=course();let face=hammer_footprint::Specification::read(
            &format!("{}\nspan,69,0.01,4",hammer_footprint::HEADER),&[c]).unwrap();
        let build=||if spatial {Bank::new_with_hammer_footprints(&[c],&board::demonstration(),192000,21600.,8,true,&face).unwrap()}
            else {Bank::new(&[c],&board::demonstration(),192000,21600.,8,true).unwrap()};
        let mut a=build();let mut b=build();
        b.configure_string_stretching(&[c],&Specification::read(&format!("{HEADER}\nlinear,69"),&[c]).unwrap()).unwrap();
        assert!(!b.has_string_stretching());
        let forces=vec![0.1;a.contact_strings.len()];
        for _ in 0..100 {a.predict();a.finish(&forces);a.commit();advance(&mut b,&forces).unwrap();b.commit();
            assert_eq!(a.q,b.q);assert_eq!(a.v,b.v);assert_eq!(a.energy(),b.energy());}
    }
}
