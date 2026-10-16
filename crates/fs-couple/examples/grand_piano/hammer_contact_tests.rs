use super::*;
use fs_material::Uniaxial;

fn residual(sites:&[Site<'_>],a:&[f64],free:&[f64],f:&[f64])->f64 {
    let n=sites.len();(0..n).map(|i| {
        let end=free[i]-(0..n).map(|j|a[i*n+j]*f[j]).sum::<f64>();
        (f[i]-sites[i].force(end).0).abs()
    }).fold(0.0_f64,f64::max)
}
fn scalar_sweeps(sites:&[Site<'_>],a:&[f64],free:&[f64],count:usize)->Vec<f64> {
    let n=sites.len();let mut f=vec![0.;n];
    for _ in 0..count {for i in 0..n {
        let s=sites[i];let offset=free[i]-(0..n).filter(|j|*j!=i).map(|j|a[i*n+j]*f[j]).sum::<f64>();
        f[i]=super::super::solve(s.law,s.history,s.start_m,offset,a[i*n+i],s.thickness_m,s.area_m2).unwrap();
    }}
    f
}

#[test]
fn simultaneous_hammer_resolves_twelve_sites_that_exhaust_scalar_sweeps() {
    // Explicit stiff/thin felt coupon, not an identified production hammer.
    // Compliance includes one shared finite hammer and twelve local responses.
    let law=WoolFelt::new(1e8,0.2,2.,2.5,0.2,0.8).unwrap();let history=law.initial_state();
    let sites=[Site {law:&law,history:&history,start_m:0.,thickness_m:1e-4,area_m2:1e-4/12.};12];
    let mut a=vec![1e-7;144];for i in 0..12 {a[12*i+i]+=1e-9;}
    let free=[2e-5;12];let scalar=scalar_sweeps(&sites,&a,&free,32);
    assert!(residual(&sites,&a,&free,&scalar)>1e-4,"fixture must expose the exhausted scalar budget");
    let mut w=Workspace::new().unwrap();let mut f=[-123.;12];
    let updates=w.solve(&sites,&a,&free,&[0.;12],&mut f).unwrap();
    assert!(updates<16 && residual(&sites,&a,&free,&f)<1e-5);
    assert!(f.iter().all(|v|(*v-f[0]).abs()<1e-8));
    assert!((f[0]-13.04413064750948).abs()<1e-7);
    // The conservative discrete work is still the exact primitive increment.
    for i in 0..12 {
        let end=free[i]-(0..12).map(|j|a[12*i+j]*f[j]).sum::<f64>();
        let energy=sites[i].area_m2*sites[i].thickness_m
            *super::super::primitive(&law,end/sites[i].thickness_m,&history);
        assert!((f[i]*end-energy).abs()<1e-12);
    }
}

#[test]
fn block_matches_converged_scalar_law_on_distinct_loading_and_unloading_histories() {
    let law=super::super::demonstration_law().unwrap();
    let histories:Vec<_>=[0.0,0.1,0.3,0.5].iter().map(|e|law.update_state(*e,&law.initial_state())).collect();
    let starts=[-0.0001,0.0004,0.0018,0.0028];
    let sites:Vec<_>=histories.iter().enumerate().map(|(i,h)|Site {
        law:&law,history:h,start_m:starts[i],thickness_m:0.008,area_m2:(i+1) as f64*1e-5,
    }).collect();
    let mut a=vec![1e-7;16];for i in 0..4 {a[5*i]+=2e-7;}
    for free in [[0.0006,0.0014,0.0025,0.001],[-0.002,-0.001,-0.001,-0.001],[0.001,0.002,0.003,0.004]] {
        let expected=scalar_sweeps(&sites,&a,&free,128);let mut f=[0.;4];
        Workspace::new().unwrap().solve(&sites,&a,&free,&[0.;4],&mut f).unwrap();
        assert!(residual(&sites,&a,&free,&f)<1e-5);
        for (x,y) in f.iter().zip(expected) {assert!((x-y).abs()<1e-5);}
        assert!(f.iter().all(|x|*x>=0.));
    }
    assert_eq!(histories[0].eps_max,0.);assert_eq!(histories[3].eps_max,0.5);
}

#[test]
fn signed_reciprocal_compliance_and_site_permutation_preserve_the_root() {
    let law=super::super::demonstration_law().unwrap();let history=law.initial_state();
    let sites=[Site {law:&law,history:&history,start_m:0.,thickness_m:0.008,area_m2:1e-4};3];
    // Positive definite, with a signed geometric cross-coupling, not abs(A).
    let a=[3e-6,-1e-6,0.5e-6,-1e-6,2e-6,-0.3e-6,0.5e-6,-0.3e-6,4e-6];
    let free=[0.001,-0.00001,0.002];let mut original=[0.;3];
    let mut w=Workspace::new().unwrap();w.solve(&sites,&a,&free,&[0.;3],&mut original).unwrap();
    let order=[2,0,1];let mut matrix=[0.;9];
    for i in 0..3 {for j in 0..3 {matrix[3*i+j]=a[3*order[i]+order[j]];}}
    let permuted=order.map(|i|sites[i]);let mut f=[0.;3];
    w.solve(&permuted,&matrix,&order.map(|i|free[i]),&[1.;3],&mut f).unwrap();
    for i in 0..3 {assert!((f[i]-original[order[i]]).abs()<1e-5);}
    assert!(residual(&sites,&a,&free,&original)<1e-5);
}

#[test]
fn block_refusal_never_publishes_and_retry_does_not_depend_on_trial_scratch() {
    let law=super::super::demonstration_law().unwrap();let h=law.initial_state();
    let site=Site {law:&law,history:&h,start_m:0.,thickness_m:0.008,area_m2:1e-4};
    let mut w=Workspace::new().unwrap();let mut out=[123.];
    for (a,free) in [(f64::NAN,0.001),(-1e-7,0.001),(1e-12,1.0)] {
        assert!(w.solve(&[site],&[a],&[free],&[0.],&mut out).is_err());assert_eq!(out,[123.]);
    }
    assert!(w.solve(&[site],&[],&[0.],&[0.],&mut out).is_err());assert_eq!(out,[123.]);
    let mut expected=[0.];
    w.solve(&[site],&[2e-7],&[0.001],&[0.],&mut out).unwrap();
    Workspace::new().unwrap().solve(&[site],&[2e-7],&[0.001],&[0.],&mut expected).unwrap();
    assert_eq!(out,expected);assert_eq!(h.eps_max,0.);
}
