use super::*;
use fs_ascent::transient::variational::intervals::IntervalTape;

fn with_domain(test: impl FnOnce(&Cx<'_>, &Domain, &BackwardEuler<'_>)) {
    let d=Domain::new().unwrap();let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:61,kernel_id:820,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let engine=BackwardEuler::uniform(&cx,&d.mesh,VolumetricHeatCapacity::declared(1000.0).unwrap()).unwrap();
        test(&cx,&d,&engine);
    });
}
#[test]
fn substep_options_have_explicit_bounds_and_preserve_default() {
    assert_eq!(options(&[]).unwrap(),(None,0.03,1));
    let args=["data.csv","--substeps","8","--model-sigma","0.02"].map(str::to_owned);
    assert_eq!(options(&args).unwrap(),(Some("data.csv"),0.02,8));
    for raw in [vec!["--substeps"],vec!["--substeps","0"],vec!["--substeps","33"],
        vec!["--substeps","1.5"],vec!["--substeps","NaN"],vec!["--substeps","2","--substeps","3"]] {
        assert!(options(&raw.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
    }
}

#[test]
fn checkpointed_conduction_matches_retained_pde_steps_and_source_derivatives() {
    with_domain(|cx,d,engine| {
        let base=ConductionWindowPolicy::new(cx,&d.mesh,&d.boundary,&[0.0,0.6],config(),&mut||false).unwrap();
        let policy=ConductionSubsteps::new(base.clone(),&[4],substep_limits(),&mut||false).unwrap();
        let n=base.dimension();let initial=vec![300.0;n];
        let family=Family {domain:d,engine,samples:Vec::new(),times:vec![0.0,0.6],dimension:n};
        let model=family.instantiate(&[1500.0,0.0],&mut||false).unwrap();
        let tape=policy.record(&model,0,0.0,0.6,&initial,&mut||false).unwrap();
        let mut current=base.expand_field(&initial).unwrap();
        for t in policy.substep_times(0).unwrap().windows(2) {
            current=engine.advance(cx,model.problem(0).unwrap(),None,&current,t[1]-t[0],config().step).unwrap().temperature;
        }
        assert_eq!(tape.endpoint(),base.gather_field(&current).unwrap());
        let seed=vec![1.0/n as f64;n];
        let got=tape.pullback(&seed,&[0.13,0.2],&mut||false).unwrap();
        assert_eq!(got.replayed_steps,8);assert_eq!(got.peak_checkpoints,3);
        let value=|source:f64| {
            let m=family.instantiate(&[source,0.0],&mut||false).unwrap();
            let t=policy.record(&m,0,0.0,0.6,&initial,&mut||false).unwrap();
            t.endpoint().iter().zip(&seed).map(|(x,b)|x*b).sum::<f64>()+0.13*source
        };
        let difference=(value(1500.1)-value(1499.9))/0.2;
        assert!((got.parameters[0]-difference).abs()<2e-7);
        let full=base.expand_field(tape.endpoint()).unwrap();
        for &(node,temperature) in d.boundary.dirichlet() {assert_eq!(full[node],temperature);}
        let expected=tape.endpoint().to_vec();
        assert!(tape.pullback(&seed,&[0.0,0.0],&mut||true).is_err());
        assert_eq!(tape.endpoint(),expected);
    });
}

#[test]
fn spatial_fit_refines_forecasts_without_adding_controls_or_defect_priors() {
    with_domain(|cx,d,engine| {
        let rows=synthetic(cx,d,engine,4).unwrap();
        let got=fit(cx,d,engine,&rows,0.03,4).unwrap();
        assert_eq!(got.report.reason,StopReason::GradNorm);
        assert_eq!(got.evaluation.controls.len(),146);
        assert_eq!(got.evaluation.window.defects.len(),108);
        assert_eq!(got.times,vec![0.0,0.25,0.5,1.0]);
        assert_eq!(got.evaluation.window.accepted_steps,12);
        assert_eq!(got.evaluation.window.replayed_steps,24);
        assert!((got.evaluation.parameters[0]-2000.0).abs()<10.0);
        assert!((got.evaluation.parameters[1]-0.08).abs()<0.01);
        assert!(got.evaluation.value<got.before*0.001);
        for field in got.fields {
            for &(node,temperature) in d.boundary.dirichlet() {assert_eq!(field[node],temperature);}
        }
    });
}

#[test]
fn spatial_forward_refinement_approaches_finer_reference_at_same_sensor_times() {
    with_domain(|cx,d,engine| {
        let coarse=synthetic(cx,d,engine,1).unwrap();
        let refined=synthetic(cx,d,engine,4).unwrap();
        let reference=synthetic(cx,d,engine,32).unwrap();
        let loss=|values:&[Reading]| values.iter().zip(&reference).map(|(a,b)| {
            assert_eq!((a.time,a.node,a.sigma),(b.time,b.node,b.sigma));(a.value-b.value).powi(2)
        }).sum::<f64>();
        assert!(loss(&coarse)>1e-5);
        assert!(loss(&refined)<0.2*loss(&coarse));
    });
}
