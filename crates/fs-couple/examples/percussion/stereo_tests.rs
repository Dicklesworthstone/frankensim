use super::*;

#[test]
fn right_microphone_preserves_playing_inputs_and_refuses_implicit_mode_changes() {
    let mut args:Vec<_>=["drum-modal-mic","128","20","--microphone-right","-0.1,0.05,0.4",
        "--strike-speed-m-s","2"].map(String::from).to_vec();
    let right=option(&mut args).unwrap();assert_eq!(right,Some([-0.1,0.05,0.4]));
    let (positional,stroke)=crate::playing::parse(args).unwrap();
    assert_eq!(positional,["drum-modal-mic","128","20"]);assert_eq!(stroke.speed_m_s,2.0);
    for family in ["splash","drum","drum-modal","drum-stretch","snare","snare-off"] {
        assert!(admit_command(right,&format!("{family}-mic")).is_ok());
        assert!(admit_command(right,family).is_err());
        assert!(admit_command(right,&format!("{family}-wav")).is_err());
    }
    assert!(admit_command(None,"drum").is_ok());
    for text in ["--microphone-right", "--microphone-right 1,2", "--microphone-right 1,2,3,4",
        "--microphone-right NaN,2,3", "--microphone-right 1,inf,3",
        "--microphone-right 1,2,3 --microphone-right 3,2,1"] {
        let mut args:Vec<_>=text.split_whitespace().map(String::from).collect();let old=args.clone();
        assert!(option(&mut args).is_err());assert_eq!(args,old);
    }
    let neck=crate::cavity::NeckOptions {radius_m:0.005,effective_length_m:0.012,
        resistance_pa_s_m3:1000.0,azimuth_rad:0.4,axial_position_m:0.08};
    assert!(crate::cavity::admit_neck_command(Some(neck),true,"snare-mic").is_err());
}

fn tetrahedron()->Boundary {
    let nodes=[[0.01,0.01,0.01],[0.01,-0.01,-0.01],[-0.01,0.01,-0.01],[-0.01,-0.01,0.01]];
    let mut triangles=Vec::new();
    for t in [[0,1,2],[0,3,1],[0,2,3],[1,3,2]] {
        let mut p=t.map(|i|nodes[i]);
        let a=[p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]];
        let b=[p[2][0]-p[0][0],p[2][1]-p[0][1],p[2][2]-p[0][2]];
        let n=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
        if n.iter().zip(p[0]).map(|(n,p)|n*p).sum::<f64>()<0.0 {p.swap(1,2);}
        triangles.push(p);
    }
    let surface=SpherePanels::from_triangles(triangles.clone()).unwrap();
    // Small closed BEM fixture with an explicitly prescribed vertical motion.
    Boundary {triangles,weights:vec![surface.normals().iter().map(|n|n[2]).collect()],state_modes:vec![0]}
}

#[test]
fn shared_bem_source_reproduces_independent_spatial_receiver_bakes() {
    let boundary=tetrahedron();
    let receivers=[Receiver::FinitePoint([0.0,0.0,0.12]),Receiver::FinitePoint([0.1,0.0,0.08])];
    let pair=bake_receivers(&boundary,&receivers).unwrap();let mut different=false;
    for (channel,&receiver) in receivers.iter().enumerate() {
        let single=bake_receivers(&boundary,&[receiver]).unwrap();
        assert_eq!(pair[channel].pressure_gain,1.0);
        assert_eq!(pair[channel].propagation_delay_s,single[0].propagation_delay_s);
        for hz in [120.0,400.0,1200.0] {
            let w=core::f64::consts::TAU*hz;
            let actual=pair[channel].filters[0].eval(w).unwrap();
            let expected=single[0].filters[0].eval(w).unwrap();
            assert!((actual-expected).abs()<1e-12*(1.0+expected.abs()));
            let other=pair[1-channel].filters[0].eval(w).unwrap();
            different|=(actual-other).abs()>1e-10;
        }
    }
    assert!(different,"distinct points must observe distinct Green responses");
    assert!(bake_receivers(&boundary,&[]).is_err());
    assert!(bake_receivers(&boundary,&[receivers[0];3]).is_err());
    assert!(bake_receivers(&boundary,&[receivers[0],Receiver::FinitePoint([0.0;3])]).is_err());
}

// Authored transfer fixtures test ONLY channel/decimator/physics plumbing. The
// independent test above exercises real BEM and the unchanged holdout gate.
fn fixture(count:usize,channel:usize)->Bake {
    let dt=1.0/f64::from(OUTPUT_RATE);
    let filters=(0..count).map(|i| {
        let gain=(i+1) as f64*0.001;
        if channel==0 {let mut f=zero_filter(dt);f.d=gain;f} else {
            DiscreteStateSpace{n:1,a:vec![0.8],b:vec![0.2],c:vec![-gain],d:0.0,e_leftover:0.0,t_s:dt}
        }
    }).collect();
    let medium=Medium::air();let delay=if channel==0 {3.0*dt}else{11.0*dt};
    Bake {filters,range_m:delay*medium.sound_speed,medium,propagation_delay_s:delay,pressure_gain:1.0}
}
fn drum(frames:usize)->Experiment {
    crate::drum((frames*SUBSTEPS) as u64,MECHANICAL_DT,true,true).unwrap()
}

#[test]
fn real_drum_advances_once_and_stereo_pcm_equals_each_mono_observer() {
    let frames=128;let mut pair=drum(frames);let mut left=drum(frames);let mut right=drum(frames);
    let count=pair.acoustics.as_ref().unwrap().state_modes.len();
    let wav=render_baked(&mut pair,frames,10.0,&[fixture(count,0),fixture(count,1)]).unwrap();
    let a=render_baked(&mut left,frames,10.0,&[fixture(count,0)]).unwrap();
    let b=render_baked(&mut right,frames,10.0,&[fixture(count,1)]).unwrap();
    assert_eq!(wav.len(),44+4*frames);assert_eq!(&wav[22..24],&2_u16.to_le_bytes());
    assert_eq!(&wav[28..32],&(OUTPUT_RATE*4).to_le_bytes());
    for (i,frame) in wav[44..].chunks_exact(4).enumerate() {
        assert_eq!(&frame[..2],&a[44+2*i..46+2*i]);assert_eq!(&frame[2..],&b[44+2*i..46+2*i]);
    }
    assert!(wav[44..].iter().any(|b|*b!=0));assert_ne!(&a[44..],&b[44..]);
    assert_eq!(pair.system.state(),left.system.state());assert_eq!(pair.system.state(),right.system.state());
    let crate::Mechanics::Prepared(system)=&pair.system else {panic!("prepared drum");};
    assert_eq!(system.samples_rendered(),(frames*SUBSTEPS) as u64);
    let crate::Mechanics::Prepared(mono)=&left.system else {panic!("prepared mono");};
    assert_eq!(system.frame(),mono.frame());
}

#[test]
fn invalid_right_receiver_or_transfer_does_not_advance_the_drum() {
    let mut experiment=drum(8);let state=experiment.system.state().to_vec();
    let receivers=[Receiver::FinitePoint([0.08,0.05,0.35]),Receiver::FinitePoint([f64::NAN,0.0,0.0])];
    assert!(render_receivers(&mut experiment,8,20.0,&receivers).is_err());
    assert_eq!(experiment.system.state(),state);
    let count=experiment.acoustics.as_ref().unwrap().state_modes.len();
    let mut bad=fixture(count,1);bad.filters[0].d=f64::NAN;
    assert!(render_baked(&mut experiment,8,20.0,&[fixture(count,0),bad]).is_err());
    assert_eq!(experiment.system.state(),state);
    let crate::Mechanics::Prepared(system)=&experiment.system else {panic!();};
    assert_eq!(system.samples_rendered(),0);
}

#[test]
fn cold_boundary_refinement_preserves_flux_source_addresses_and_original_geometry() {
    let mut b=tetrahedron();b.weights.push(vec![1.0,2.0,-3.0,4.0]);b.state_modes.push(17);
    let triangles=b.triangles.clone();let weights=b.weights.clone();let gate=CancelGate::new_clock_free();
    let original=SpherePanels::from_triangles(triangles.clone()).unwrap();
    let flux=|rows:&[Vec<f64>],areas:&[f64]|rows.iter().map(|r|r.iter().zip(areas).map(|(b,a)|b*a).sum::<f64>()).collect::<Vec<_>>();
    let expected=flux(&b.weights,original.areas());
    for levels in 0..=2 {
        let refined=curved_aperture::uniform_refinement(&b,levels,64,&gate).unwrap();
        assert_eq!(refined.triangles.len(),4*4_usize.pow(levels));assert_eq!(refined.state_modes,b.state_modes);
        let surface=SpherePanels::from_triangles(refined.triangles.clone()).unwrap();
        for (actual,want) in flux(&refined.weights,surface.areas()).iter().zip(&expected) {assert!((actual-want).abs()<1e-16);}
        for p in triangles.iter().flatten() {assert!(refined.triangles.iter().flatten().any(|q|q==p));}
        for (i,normal) in surface.normals().iter().enumerate() {
            let parent=i/(4_usize.pow(levels));
            assert!(normal.iter().zip(&original.normals()[parent]).map(|(a,b)|a*b).sum::<f64>()>1.0-1e-12);
        }
    }
    assert_eq!(b.triangles,triangles);assert_eq!(b.weights,weights);
    assert!(curved_aperture::uniform_refinement(&b,2,63,&gate).is_err());
    gate.request();assert!(curved_aperture::uniform_refinement(&b,1,64,&gate).is_err());
}

#[test]
fn wider_radiation_band_uses_real_bem_and_independently_audited_receiver_filters() {
    let mut b=tetrahedron();
    // Prescribed breathing of a small closed polyhedron, not an eigensolve or a
    // measured cymbal. It tests the full wider-band BEM -> digital filter path.
    b.weights[0].fill(1.0);
    let spec=radiation_spec::Spec {band_hz:[40.0,8000.0],training_intervals:64,max_order:24,
        subdivisions:2,max_panels:64,max_dense_work:1_000_000_000};
    let receivers=[Receiver::FinitePoint([0.0,0.0,0.12]),Receiver::FinitePoint([0.1,0.0,0.08])];
    let gate=CancelGate::new_clock_free();let result=bake_receivers_with_spec(&b,&receivers,spec,&gate).unwrap();
    assert_eq!(result.len(),2);assert_eq!(result[0].filters.len(),1);
    let refined=curved_aperture::uniform_refinement(&b,2,64,&gate).unwrap();
    let surface=SpherePanels::from_triangles(refined.triangles.clone()).unwrap();
    let radius=b.triangles.iter().flatten().map(|p|p.iter().map(|x|x*x).sum::<f64>().sqrt()).fold(0.0_f64,f64::max);
    let medium=Medium::air();let mut compared=0;
    // New frequencies are absent from training, order selection AND audit.
    for hz in [2207.0,3911.0,6123.0,7777.0] {
        let w=core::f64::consts::TAU*hz;let fields=acceleration_fields(&refined.weights,w);
        let solutions=solve_radiation_batch(&surface,w/medium.sound_speed,medium,&[&fields[0]],Formulation::BurtonMiller).unwrap();
        for (channel,&receiver) in receivers.iter().enumerate() {
            let truth=receiver_response(&surface,&solutions[0],medium,receiver,radius,result[channel].propagation_delay_s).unwrap();
            let fitted=result[channel].filters[0].eval(w).unwrap().conj();
            assert!((fitted-truth).abs()<0.15*truth.abs(),"out-of-grid frequency {hz}");compared+=1;
        }
    }
    assert_eq!(compared,8);
}

#[test]
fn acoustic_budget_cancellation_and_neck_band_refusals_never_advance_mechanics() {
    let mut e=drum(4);let before=e.system.state().to_vec();let receiver=[Receiver::FinitePoint([0.08,0.05,0.35])];
    let gate=CancelGate::new_clock_free();
    let tiny=radiation_spec::Spec {max_dense_work:1,..radiation_spec::Spec::default()};
    assert!(render_receivers_with_spec(&mut e,4,20.0,&receiver,tiny,&gate).is_err());
    assert_eq!(e.system.state(),before);
    gate.request();assert!(render_receivers_with_spec(&mut e,4,20.0,&receiver,radiation_spec::Spec::default(),&gate).is_err());
    assert_eq!(e.system.state(),before);
    let neck=crate::cavity::NeckOptions {radius_m:0.003,effective_length_m:0.008,
        resistance_pa_s_m3:1000.0,azimuth_rad:0.4,axial_position_m:0.08};
    let e=crate::drum_with_radiation(16,MECHANICAL_DT,false,true,None,false,crate::Stroke::default(),
        true,Some(neck),None,None,&[],0.0,None,false).unwrap();
    let broad=radiation_spec::Spec {band_hz:[40.0,8000.0],..radiation_spec::Spec::default()};
    assert!(broad.admit_necks(&e).is_err());radiation_spec::Spec::default().admit_necks(&e).unwrap();
}
