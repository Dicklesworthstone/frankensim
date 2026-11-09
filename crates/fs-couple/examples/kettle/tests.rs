use super::{model::{self,Model},spec::Spec};
use fs_exec::CancelGate;
fn spec() -> Spec { Spec::read(include_str!("illustrative-parabolic.fsk")).unwrap() }
fn uniform() -> Spec {
    let mut s=spec();
    s.head.radial_intervals=2; s.head.azimuths=16;
    s.head_band_hz=[0.0,250.0]; s.air_maximum_hz=1.0;
    s.air_mesh=[2,0]; s.creep.clear();
    s
}
#[test]
fn specification_preserves_geometry_and_refuses_incomplete_or_duplicate_data() {
    let s=spec();
    assert_eq!(s.head.radius_m,0.3302); assert_eq!(s.profile.len(),9);
    let text=include_str!("illustrative-parabolic.fsk");
    for changed in [text.replace("head_mesh,4,32",""),format!("{text}\nrate_hz,768000"),
        text.replace("station,0,0","station,0,0.01"), text.replace("air,1.2,343,0","air,1.2,NaN,0"),
        text.replace("rate_hz,768000","rate_hz,0"),format!("{text}\nunknown,2")] {
        assert!(Spec::read(&changed).is_err());
    }
}
#[test]
fn complete_mallet_disk_cannot_cross_polygon_between_sampling_points() {
    let s=uniform(); let gate=CancelGate::new_clock_free();
    let m=Model::new(&s,1,&gate).unwrap();
    assert_eq!(model::footprint(&m.film,[0.0,0.0],0.01).unwrap().len(),4);
    assert!(model::footprint(&m.film,[s.head.radius_m-0.001,0.0],0.01).is_err());
}
#[test]
fn uniform_overlap_is_signed_actual_head_area_and_one_bulk_storage() {
    let s=uniform(); let gate=CancelGate::new_clock_free();
    let m=Model::new(&s,1,&gate).unwrap();
    let (air,overlap)=model::interface(&m.film,&m.modes,&m.air,&gate).unwrap();
    assert_eq!(air.omegas,vec![0.0]);
    assert_eq!(m.cavity.total_modes(),1+m.modes.len()); // no fake zero-mode inertia
    assert_eq!(overlap[0],0.0); // mallet never directly pumps the air
    for (i,mode) in m.modes.iter().enumerate() {
        let area=m.film.modal_area(&mode.phi).unwrap();
        assert!((overlap[i+1]+area).abs()<1e-12);
    }
    let mut state=m.mechanics.state().to_vec(); state[2]=1e-5;
    let expected=-s.air[0]*s.air[1].powi(2)*overlap[1]*state[2]/m.air.volume_m3();
    let mut p=[0.0]; m.cavity.pressures_into(&state,&mut p).unwrap();
    assert!((p[0]-expected).abs()<1e-10*expected.abs().max(1.0));
}
#[test]
fn nonlinear_head_air_and_felt_preserve_history_on_cancel_and_force_refusal() {
    let s=uniform(); let gate=CancelGate::new_clock_free();
    let mut a=Model::new(&s,800,&gate).unwrap(); let mut b=Model::new(&s,800,&gate).unwrap();
    let initial=a.mechanics.stored_energy_j(); let mut loss=0.0; let mut moved=false;
    for tick in 0..800 {
        if tick==200 {
            let before=a.mechanics.state().to_vec(); let history=a.mechanics.felt_history(0).unwrap();
            let cancelled=CancelGate::new_clock_free(); cancelled.request();
            assert!(a.mechanics.step(&a.force,&cancelled).is_err());
            let mut bad=a.force.clone(); bad[0]=f64::NAN;
            assert!(a.mechanics.step(&bad,&gate).is_err());
            assert_eq!(a.mechanics.state(),before);
            assert_eq!(a.mechanics.felt_history(0).unwrap(),history);
        }
        let frame=a.mechanics.step(&a.force,&gate).unwrap();
        b.mechanics.step(&b.force,&gate).unwrap();
        assert_eq!(a.mechanics.state(),b.mechanics.state());
        loss+=frame.dissipated_energy_j;
        assert!((a.mechanics.stored_energy_j()+loss-initial).abs()<1e-6);
        moved|=a.head_displacement().abs()>1e-12;
    }
    assert!(moved); assert!(a.mechanics.felt_history(0).unwrap().eps_max>0.0);
    let before=a.mechanics.state().to_vec();
    assert!(a.mechanics.step(&a.force,&gate).is_err()); assert_eq!(a.mechanics.state(),before);
}
#[test]
fn hysteresis_free_felt_law_strikes_the_head() {
    // F = K x^p unloading on its own curve (crush 0, q = p), as Rhaouti, Chaigne and Joly's mallet.
    let mut s=uniform(); s.felt=[100000.0,0.2,2.2,2.2,0.0,0.7,0.0];
    let gate=CancelGate::new_clock_free();
    let mut m=Model::new(&s,800,&gate).unwrap();
    let initial=m.mechanics.stored_energy_j(); let mut loss=0.0;
    for _ in 0..800 { loss+=m.mechanics.step(&m.force,&gate).unwrap().dissipated_energy_j; }
    assert!(m.mechanics.felt_history(0).unwrap().eps_max>0.0);
    assert!((m.mechanics.stored_energy_j()+loss-initial).abs()<1e-6);
    s.felt[3]=3.0; // crush 0 is the law's limit only when q = p
    assert!(Model::new(&s,800,&gate).is_err());
}
#[test]
fn an_explicit_zero_exterior_load_preserves_the_unloaded_model_exactly() {
    use fs_couple::render::plate::impact::radiation::{Model as Load,Pole};
    let s=uniform();let gate=CancelGate::new_clock_free();
    let mut a=Model::new(&s,40,&gate).unwrap();
    let mut b=Model::with_exterior(&s,40,&gate,|_,modes|Ok(Some(Load {
        ports:modes.len(),poles:vec![Pole {omega:1000.0,zeta:0.2,coupling:vec![0.0;modes.len()]}]
    }))).unwrap();
    for _ in 0..40 {
        a.mechanics.step(&a.force,&gate).unwrap(); b.mechanics.step(&b.force,&gate).unwrap();
        assert_eq!(a.mechanics.state(),b.mechanics.state());
    }
}
#[test]
fn passive_exterior_memory_joins_the_same_impact_and_rollback() {
    use fs_couple::render::plate::impact::radiation::{Model as Load,Pole};
    let s=uniform();let gate=CancelGate::new_clock_free();
    // Declared passive fixture tests composition; it is not a BEM sound model.
    let mut a=Model::with_exterior(&s,800,&gate,|_,modes|Ok(Some(Load {
        ports:modes.len(),poles:vec![Pole {omega:1000.0,zeta:0.2,coupling:vec![100.0;modes.len()]}]
    }))).unwrap();
    let initial=a.mechanics.stored_energy_j();let mut loss=0.0;
    for _ in 0..800 { loss+=a.mechanics.step(&a.force,&gate).unwrap().dissipated_energy_j; }
    assert!((a.mechanics.stored_energy_j()+loss-initial).abs()<1e-6);
    assert!(a.mechanics.radiation_observation().unwrap().stored_energy_j>0.0);
    let state=a.mechanics.state().to_vec();let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(a.mechanics.step(&a.force,&cancelled).is_err());assert_eq!(a.mechanics.state(),state);
}
fn rim_nodes(film:&fs_plate::shell::head::TensionedDisk)->Vec<usize> {
    let mut fixed:Vec<_>=model::rim(film).iter().flat_map(|&(a,b)|[a,b]).collect();
    fixed.sort_unstable(); fixed.dedup(); fixed
}
fn pet_pack()->super::film_loss::FilmLoss {
    super::film_loss::FilmLoss::load(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk")).unwrap()
}
#[test]
fn the_film_loss_acts_on_the_bending_share_of_the_timpani_head() {
    // The timpani card's simply supported head (w = 0, free slope): each mode's bending share phi^T K_b phi / lambda
    // against the simply supported plate's own: D k^2 (1 - 2 (1 - nu) / (k a)^2) / (T + D k^2), k^2 = rho h omega^2 / T
    // at the mode's in-vacuo frequency. The (1 - 2 (1 - nu) / (k a)^2) is the curved free-slope edge's Gaussian-
    // curvature energy, (1 - nu) / a times the edge integral of (dw/dn)^2, taken off D times the integral of
    // (lap w)^2. On the same head at 24 radial intervals (48 azimuths, as the card) every retained mode agrees
    // within 3 % and (0,1) and the (1,1) pair within 2 %; the card's own 6 radial intervals (u-head-mesh) are coarse
    // for curvature, so its shares are bounded and printed, and they are the ones the model's loss uses.
    use fs_plate::shell::{head::TensionedDisk, profile::ProfileBudget};
    let s=Spec::load(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../data/instruments/descriptions/timpani/timpani-2b.fsk")).unwrap();
    let ratios=|radial:usize|->Vec<f64> {
        let mut head=s.head; head.radial_intervals=radial;
        let film=TensionedDisk::new(head,ProfileBudget {max_nodes:2048,max_triangles:4096,max_feature_evaluations:100_000}).unwrap();
        let modes=fs_modal::slice_window(&film.model.k,&film.model.m,
            ((std::f64::consts::TAU*s.head_band_hz[0]).powi(2),(std::f64::consts::TAU*s.head_band_hz[1]).powi(2)),
            &fs_plate::SliceOptions::default()).unwrap().modes;
        let shares=super::film_loss::bending_shares(&film,&modes,&rim_nodes(&film)).unwrap();
        let d=head.young_pa*head.thickness_m.powi(3)/(12.0*(1.0-head.poisson*head.poisson));
        modes.iter().zip(&shares).enumerate().map(|(n,(mode,share))| {
            let k2=head.density_kg_m3*head.thickness_m*mode.lambda/head.tension_n_m;
            let plate=d*k2*(1.0-2.0*(1.0-head.poisson)/(k2*head.radius_m*head.radius_m))/(head.tension_n_m+d*k2);
            eprintln!("radial {radial}: mode {n}: {:.3} Hz, bending share {share:.6e}, simply supported plate {plate:.6e}, ratio {:.4}",
                mode.lambda.sqrt()/std::f64::consts::TAU,share/plate);
            assert!(*share>0.0 && *share<2e-3,"radial {radial} mode {n}: share {share}");
            share/plate
        }).collect()
    };
    let card=ratios(s.head.radial_intervals);
    assert_eq!((s.head.radial_intervals,s.head.azimuths,card.len()),(6,48,17));
    let fine=ratios(24);
    assert_eq!(fine.len(),17);
    for (n,ratio) in fine.iter().enumerate() {
        assert!((ratio-1.0).abs()<if n<3 {0.02} else {0.03},"mode {n}: ratio {ratio} at 24 radial intervals");
    }
}
#[test]
fn without_the_flag_the_head_keeps_the_cards_ratio_and_with_it_the_loss_is_bending_only() {
    let gate=CancelGate::new_clock_free();
    let mut s=uniform();
    assert!(s.head_damping_ratio>0.0);
    let m=Model::new(&s,1,&gate).unwrap();
    let card:Vec<f64>=m.modes.iter().map(|mode|2.0*s.head_damping_ratio*mode.lambda.sqrt()).collect();
    assert_eq!(m.head_damping_per_s,card);
    s.head_damping_ratio=0.0;
    let m=Model::new(&s,1,&gate).unwrap();
    assert!(m.head_damping_per_s.iter().all(|d|*d==0.0));
    let pet=pet_pack();
    s.head_film_loss=Some(std::sync::Arc::new(pet_pack()));
    let m=Model::new(&s,1,&gate).unwrap();
    let shares=super::film_loss::bending_shares(&m.film,&m.modes,&rim_nodes(&m.film)).unwrap();
    assert_eq!(m.head_damping_per_s.len(),m.modes.len());
    for ((rate,mode),share) in m.head_damping_per_s.iter().zip(&m.modes).zip(&shares) {
        let omega=mode.lambda.sqrt();
        assert_eq!(*rate,pet.tan_delta(omega/std::f64::consts::TAU).unwrap()*share*omega);
        // Far below the whole stiffness's tan delta * omega: the prestress carries none of it.
        assert!(*rate>0.0 && *rate<1e-2*pet.tan_delta(omega/std::f64::consts::TAU).unwrap()*omega);
    }
    s.head_damping_ratio=0.001;
    assert_eq!(Model::new(&s,1,&gate).err().unwrap().to_string(),super::film_loss::VISCOUS_RATIO_REFUSAL);
}
#[test]
fn film_loss_receipts_and_flags_are_explicit() {
    let gate=CancelGate::new_clock_free();
    let mut s=uniform(); s.head_damping_ratio=0.0;
    let m=Model::new(&s,1,&gate).unwrap();
    let (rates,lines)=pet_pack().damping(&m.film,&m.modes,&rim_nodes(&m.film)).unwrap();
    assert_eq!((rates.len(),lines.len()),(m.modes.len(),m.modes.len()));
    for (n,line) in lines.iter().enumerate() {
        for part in [format!("mode {n} in vacuo"),"bending share".into(),"tan delta".into(),"293.15 K".into(),
            "on bending only (prestress carries no loss)".into(),"head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk".into(),
            "measured on another film than this head".into()] {
            assert!(line.contains(&part),"{line} lacks {part}");
        }
    }
    let pack=concat!(env!("CARGO_MANIFEST_DIR"),"/../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
    let mut args:Vec<String>=["a.fsk","1","--head-film-loss",pack,"--wav-float32"].map(String::from).to_vec();
    assert!(super::film_loss::option(&mut args).unwrap().is_some());
    assert!(super::audio::float32_option(&mut args).unwrap());
    assert_eq!(args,["a.fsk","1"]);
    assert!(super::film_loss::option(&mut args).unwrap().is_none());
    assert!(!super::audio::float32_option(&mut args).unwrap());
    let mut missing:Vec<String>=["a.fsk","1","--head-film-loss"].map(String::from).to_vec();
    assert!(super::film_loss::option(&mut missing).is_err());
    let mut twice:Vec<String>=["--wav-float32","--wav-float32"].map(String::from).to_vec();
    assert!(super::audio::float32_option(&mut twice).is_err());
    assert!(pet_pack().tan_delta(10.0).is_err()); // below the curve's 20 Hz band
}
