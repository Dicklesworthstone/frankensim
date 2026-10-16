use super::*;
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_exec::{Budget,CancelGate,ExecMode,StreamKey};

fn context<T>(gate:&CancelGate,f:impl FnOnce(&Cx<'_>)->T)->T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena|f(&Cx::new(gate,arena,
        StreamKey {seed:41,kernel_id:721,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn config()->ViewFactorConfig {
    ViewFactorConfig {rays_per_surface:4096,max_triangle_tests:2_000_000,
        max_balance_iterations:512,max_sampling_change:0.04,max_factor_adjustment:0.04}
}
// Test-only independent low-discrepancy stream; production callers own theirs.
fn sample(surface:usize,ray:u32)->[f64;5] {
    [2_u32,3,5,7,11].map(|base| {
        let mut index=ray+1+surface as u32*104729;
        let mut factor=1.0;let mut value=0.0;
        while index>0 {factor/=f64::from(base);value+=f64::from(index%base)*factor;index/=base;}
        value
    })
}
fn append_box(positions:&mut Vec<[f64;3]>,triangles:&mut Vec<[u32;3]>,owners:&mut Vec<Option<usize>>,
    lo:[f64;3],hi:[f64;3],inward:bool,obstacle:bool) {
    let offset=positions.len() as u32;
    let [a,b,c]=lo;let [x,y,z]=hi;
    positions.extend([[a,b,c],[x,b,c],[x,y,c],[a,y,c],[a,b,z],[x,b,z],[x,y,z],[a,y,z]]);
    for (owner,q) in [[0,4,7,3],[1,2,6,5],[0,1,5,4],[3,7,6,2],[0,3,2,1],[4,5,6,7]].iter().enumerate() {
        for mut t in [[q[0],q[1],q[2]],[q[0],q[2],q[3]]] {
            if inward {t.swap(1,2);}
            triangles.push(t.map(|i|i+offset));owners.push(Some(if obstacle{6}else{owner}));
        }
    }
}
fn enclosure()->(Vec<[f64;3]>,Vec<[u32;3]>,Vec<Option<usize>>) {
    let (mut p,mut t,mut o)=(Vec::new(),Vec::new(),Vec::new());
    append_box(&mut p,&mut t,&mut o,[0.0;3],[2.0,1.0,1.0],true,false);(p,t,o)
}
// Independent closed-form integral for opposed, aligned equal rectangles.
fn opposed(w:f64,h:f64,d:f64)->f64 {
    let x=w/d;let y=h/d;
    2.0/(core::f64::consts::PI*x*y)*(
        0.5*((1.0+x*x)*(1.0+y*y)/(1.0+x*x+y*y)).ln()
        +x*(1.0+y*y).sqrt()*(x/(1.0+y*y).sqrt()).atan()
        +y*(1.0+x*x).sqrt()*(y/(1.0+x*x).sqrt()).atan()
        -x*x.atan()-y*y.atan())
}
#[test]
fn rectangular_enclosure_matches_independent_integral_and_closes() {
    let (p,t,o)=enclosure();
    context(&CancelGate::new_clock_free(),|cx| {
        let r=estimate_view_factors(cx,&p,&t,&o,6,config(),sample).unwrap();
        for (i,j,w,h,d) in [(0,1,1.0,1.0,2.0),(2,3,2.0,1.0,1.0),(4,5,2.0,1.0,1.0)] {
            assert!((r.raw_factors()[i][j]-opposed(w,h,d)).abs()<0.02);
            assert!((r.factors()[i][j]-opposed(w,h,d)).abs()<0.02);
        }
        for i in 0..6 {
            assert_eq!(r.counts()[i].iter().sum::<u32>(),config().rays_per_surface);
            assert_eq!(r.factors()[i][i],0.0);
            assert!((r.factors()[i].iter().sum::<f64>()-1.0).abs()<2e-13);
            for j in 0..6 {
                assert!((r.areas()[i]*r.factors()[i][j]-r.areas()[j]*r.factors()[j][i]).abs()<1e-13);
            }
        }
        assert_eq!(r.triangle_tests(),6*4096*11);
        assert!(r.adjustment()>0.0 && r.adjustment()<config().max_factor_adjustment);
    });
}
#[test]
fn nearest_opaque_obstacle_blocks_exchange_and_unassigned_hits_refuse() {
    let (mut p,mut t,mut o)=enclosure();
    context(&CancelGate::new_clock_free(),|cx| {
        let clear=estimate_view_factors(cx,&p,&t,&o,6,config(),sample).unwrap();
        append_box(&mut p,&mut t,&mut o,[0.8,0.3,0.3],[1.2,0.7,0.7],false,true);
        let shielded=estimate_view_factors(cx,&p,&t,&o,7,config(),sample).unwrap();
        assert!(shielded.raw_factors()[0][6]>0.03);
        assert!(shielded.raw_factors()[0][1]<clear.raw_factors()[0][1]);
        for owner in &mut o[12..] {*owner=None;}
        assert!(matches!(estimate_view_factors(cx,&p,&t,&o,6,config(),sample),Err(ViewFactorError::UnassignedHit{..})));
    });
}
#[test]
fn replay_scale_and_geometry_binding() {
    let (p,t,o)=enclosure();
    context(&CancelGate::new_clock_free(),|cx| {
        let a=estimate_view_factors(cx,&p,&t,&o,6,config(),sample).unwrap();
        let b=estimate_view_factors(cx,&p,&t,&o,6,config(),sample).unwrap();
        assert_eq!(a.factors(),b.factors());assert_eq!(a.geometry_identity(),b.geometry_identity());
        let scaled:Vec<_>=p.iter().map(|p|p.map(|v|v*32.0)).collect();
        let b=estimate_view_factors(cx,&scaled,&t,&o,6,config(),sample).unwrap();
        assert_eq!(a.counts(),b.counts());assert_eq!(a.factors(),b.factors());
        assert_ne!(a.geometry_identity(),b.geometry_identity());
        for (x,y) in a.areas().iter().zip(b.areas()) {assert_eq!(*x*1024.0,*y);}
    });
}
#[test]
fn escaping_rays_and_reversed_faces_never_become_ambient_or_self_factors() {
    let (p,mut t,mut o)=enclosure();
    t.truncate(10);o.truncate(10);
    // Keep every patch nonempty while removing the ceiling.
    for owner in &mut o {if *owner==Some(4){*owner=Some(0);}}
    o[0]=Some(4);o[1]=Some(4);
    context(&CancelGate::new_clock_free(),|cx| {
        assert!(matches!(estimate_view_factors(cx,&p,&t,&o,5,config(),sample),Err(ViewFactorError::Escape{..})));
        let (p,mut t,o)=enclosure();t[0].swap(1,2);
        assert!(estimate_view_factors(cx,&p,&t,&o,6,config(),sample).is_err());
    });
}
#[test]
fn work_cancel_samples_and_projection_policies_fail_closed() {
    let (p,t,o)=enclosure();let gate=CancelGate::new_clock_free();gate.request();
    context(&gate,|cx|assert_eq!(estimate_view_factors(cx,&p,&t,&o,6,config(),sample).unwrap_err(),ViewFactorError::Interrupted));
    context(&CancelGate::new_clock_free(),|cx| {
        let mut c=config();c.max_triangle_tests=1;
        let calls=std::cell::Cell::new(0);
        assert!(matches!(estimate_view_factors(cx,&p,&t,&o,6,c,|_,_|{calls.set(calls.get()+1);[0.5;5]}),Err(ViewFactorError::Budget{..})));
        assert_eq!(calls.get(),0);
        assert!(matches!(estimate_view_factors(cx,&p,&t,&o,6,config(),|_,_|[f64::NAN;5]),Err(ViewFactorError::Invalid(_))));
        let mut c=config();c.max_factor_adjustment=1e-15;
        assert!(matches!(estimate_view_factors(cx,&p,&t,&o,6,c,sample),Err(ViewFactorError::Adjustment{..})));
        assert!(matches!(balance(cx,&[1.0,2.0],&[vec![0.0,1.0],vec![1.0,0.0]],8),Err(ViewFactorError::Balance{..})));
    });
}
#[test]
fn cancellation_from_the_sampler_is_observed_before_another_ray_can_publish() {
    let (p,t,o)=enclosure();let gate=CancelGate::new_clock_free();
    context(&gate,|cx| {
        let result=estimate_view_factors(cx,&p,&t,&o,6,config(),|s,r|{if r==3{gate.request();}sample(s,r)});
        assert_eq!(result.unwrap_err(),ViewFactorError::Interrupted);
    });
}
