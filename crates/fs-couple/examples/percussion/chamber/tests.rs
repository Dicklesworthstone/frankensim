use super::*;
const INPUT:&str=include_str!("../estimated-waisted-chamber.fscp");
#[test]
fn profile_input_is_bounded_explicit_and_preserves_other_playing_arguments() {
    let profile=Spec::read(INPUT).unwrap();let drum=drum_spec::Spec::reference();
    profile.admit(&drum,2e-6,false,false).unwrap();
    assert!(profile.volume_m3()<drum.volume_m3());
    for invalid in [INPUT.replace("station,0,", "station,0.01,"),
        INPUT.replace("0.12,0.1275", "0.12,0.1"),
        INPUT.replace("mesh,6,2", "mesh,6,2\nmesh,6,2"),
        INPUT.replace("band_hz,1100", "band_hz,NaN"),
        INPUT.replace("1000000,8", "1000000,9"),
        INPUT.replace("0.08255", "0.041275")] {
        assert!(Spec::read(&invalid).is_err());
    }
    let mut args=vec!["drum".into(),"--chamber-profile".into(),"input.fscp".into(),"--strike-speed-m-s".into(),"4".into()];
    assert_eq!(option(&mut args).unwrap().as_deref(),Some("input.fscp"));
    assert_eq!(args,["drum","--strike-speed-m-s","4"]);
    for input in [vec!["--chamber-profile"],vec!["--chamber-profile","--cavity-modes"],
        vec!["--chamber-profile","a","--chamber-profile","b"]] {
        let mut args:Vec<String>=input.into_iter().map(String::from).collect();let before=args.clone();
        assert!(option(&mut args).is_err());assert_eq!(args,before);
    }
    assert!(profile.admit(&drum,2e-6,true,false).is_err());
    assert!(profile.admit(&drum,2e-6,false,true).is_err());
    assert!(profile.admit(&drum,0.01,false,false).is_err());
    assert!(profile.admit(&drum_spec::Spec {depth_m:0.17,..drum},2e-6,false,false).is_err());
}
#[test]
fn profile_pressure_uses_actual_volume_and_the_unchanged_signed_head_overlap() {
    let drum=drum_spec::Spec {radial_intervals:2,azimuths:16,..drum_spec::Spec::reference()};
    let (films,modes)=drum.prepare(2e-6,false).unwrap();
    let profile=Spec::read(INPUT).unwrap();let gate=CancelGate::new_clock_free();
    let air=profile.air(AcousticMedium {rho0:1.2,c0:343.0},&gate).unwrap();let count=air.modes().len();
    let rows=head_overlaps(&air,&films,&modes,1+modes.iter().map(Vec::len).sum::<usize>(),&gate).unwrap();
    assert!((air.volume_m3()/profile.volume_m3()-1.0).abs()<1e-14);
    assert_eq!(air.modes()[0].norm_m3,air.volume_m3());
    assert!(rows[..count].iter().all(|c|*c==0.0));
    let mut at=1;
    for (head,(film,pairs)) in films.iter().zip(&modes).enumerate() {
        for pair in pairs {
            let area=film.modal_area(&pair.phi).unwrap();
            assert!((rows[at*count]-if head==0 {-area}else{area}).abs()<1e-12);at+=1;
        }
    }
    assert!(air.values_at([0.15,0.0,0.08255]).is_err(),"inside old cylinder but outside supplied waist");
}
fn experiment(prepared:bool, profile:Option<&Spec>) -> super::super::Experiment {
    let drum=drum_spec::Spec {radial_intervals:2,azimuths:16,..drum_spec::Spec::reference()};
    super::super::drum_with_chamber(128,2e-6,false,prepared,None,false,
        super::super::Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])},true,None,
        Some(drum),None,&[],0.0,None,false,None,&super::super::mallets::Selection::default(),
        &super::super::shaft_playing::Selection::default(),None,None,profile).unwrap()
}
#[test]
fn both_existing_time_images_accept_the_profile_and_retain_the_inertial_strike() {
    let profile=Spec::read(INPUT).unwrap();
    for prepared in [false,true] {
        let mut e=experiment(prepared,Some(&profile));let gate=CancelGate::new_clock_free();
        let initial=0.5*e.system.state().chunks_exact(2).map(|q|q[1]*q[1]).sum::<f64>();
        let mut energy=initial;
        for _ in 0..128 {
            let f=e.system.step(&e.force,&gate).unwrap();
            assert!(f.balance_residual_j.abs()<1e-8);
            assert_eq!(f.supplied_work_j,0.0);energy=f.stored_energy_j;
        }
        assert!(e.system.state()[2..].iter().any(|v|v.abs()>1e-12));
        assert!(energy.is_finite());
        assert!(energy<=initial+1e-8);
        assert_eq!(e.pressure.as_ref().unwrap().volume_m3,profile.volume_m3());
        assert!(e.air.as_ref().unwrap().points(e.system.state()).unwrap().0.is_finite());
    }
}
