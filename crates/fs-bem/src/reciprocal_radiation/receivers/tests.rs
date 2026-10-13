use super::*;
use crate::near_field::FirstOrder;

fn budget()->Budget {Budget {maximum_unknowns:128,maximum_cache_bytes:4*1024*1024,
    maximum_work:100_000_000,angular_refinement:0}}
fn medium()->Medium {Medium {density:1.2,sound_speed:343.}}
fn options()->Options {Options {relative_tolerance:1e-8,maximum_depth:14,maximum_kernel_evaluations:2_000_000}}
fn surface()->SpherePanels {
    let p=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    SpherePanels::from_triangles([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]
        .map(|t|t.map(|i|p[i])).to_vec()).unwrap()
}
fn manufactured<'p,'s>(owner:&'p Prepared<'s>)->Solution<'p,'s> {
    let zero=vec![C64::ZERO;owner.facets.len()];
    let mut solution=owner.solve_batch(0.7,medium(),&[&zero],||false).unwrap().remove(0);
    // Observer-only manufactured layer densities. Overwrite the coefficients
    // intentionally; the zero solution's panel means/energy are NOT evidence
    // for this field. The observer must read the P1/P0 densities, not those means.
    for f in &owner.facets {for (id,p) in f.ids.iter().zip(f.triangle) {
        solution.coefficients[*id]=C64::new(1.+0.7*p[0]-0.2*p[1]+0.5*p[2],
            -0.3+0.4*p[0]+0.1*p[1]-0.6*p[2]);
    }}
    for (i,v) in [C64::new(0.2,0.1),C64::new(-0.3,0.4),C64::new(0.5,-0.2),C64::new(0.1,0.3)]
        .into_iter().enumerate() {solution.coefficients[owner.vertices+i]=v;}
    solution
}

#[test]
fn g1_close_linear_field_matches_independent_potential_and_velocity_integrals() {
    let mesh=surface();let owner=Prepared::new(&mesh,0.7,budget(),||false).unwrap();
    let solution=manufactured(&owner);let positions=[[0.22,0.18,-0.08]];
    let geometry=owner.receivers(&positions,0.01).unwrap();
    assert!((geometry.clearances_m()[0]-0.08).abs()<1e-14);
    let rows=geometry.prepare(0.7,medium(),options(),true,||false).unwrap();
    let actual=rows.evaluate_velocity(&solution).unwrap();
    // Independent 180x180 tensor-Gauss integral over the original tetrahedron,
    // checked against 120x120 and a spatial finite difference of scalar pressure.
    // No boundary-solve or measured microphone claim belongs to this fixture.
    let expected=C64::new(-0.033093146205789736,0.12137315152119399);
    let velocity=[C64::new(0.0004893232084501287,-0.00013259287357377374),
        C64::new(0.00013114262310612862,-0.00007735705949176258),
        C64::new(0.0009897789536010913,0.0003204154991354306)];
    assert!((actual.scalar.pressure[0]-expected).abs()<1e-9);
    for (a,b) in actual.particle_velocity_m_s[0].iter().zip(velocity) {assert!((*a-b).abs()<1e-10);}
    assert!(rows.kernel_evaluations()>4*80,"near geometry must refine");
    assert!(actual.scalar.quadrature_error_estimate_pa[0]<1e-6);
    assert!(actual.quadrature_error_estimate_m_s[0].iter().all(|e|*e<1e-6));
    let mean=rows.rows.evaluate(|i| {
        let f=&owner.facets[i];let p=f.ids.iter().fold(C64::ZERO,|s,&a|s+solution.coefficients[a]).scale(1./3.);
        ([p;3],solution.coefficients[owner.vertices+i])
    },false).unwrap().scalar.pressure[0];
    assert!((mean-expected).abs()>0.04,"flattening P1 pressure is a materially wrong receiver");
    // Opposite ideal cardioids sum to the colocated pressure; their reactive
    // response comes from this field velocity, not an assumed source direction.
    let a=FirstOrder::new(0.5,[0.,0.,1.]).unwrap();let b=FirstOrder::new(0.5,[0.,0.,-1.]).unwrap();
    let front=a.observe(expected,velocity,medium()).unwrap();let back=b.observe(expected,velocity,medium()).unwrap();
    assert!((front+back-expected).abs()<1e-14);assert!((front-back).abs()>0.1);
}

#[test]
fn g3_real_source_gradient_matches_pressure_differences_and_reuses_each_receiver_row() {
    let mesh=SpherePanels::icosphere(1.,0).unwrap();
    let owner=Prepared::new(&mesh,0.7,budget(),||false).unwrap();
    let v:Vec<_>=mesh.normals().iter().map(|n|C64::new(n[0]+0.2,n[2]-0.4*n[1])).collect();
    let scaled:Vec<_>=v.iter().map(|v|*v*C64::new(-0.4,0.3)).collect();
    let sources=owner.solve_batch(0.5,medium(),&[&v,&scaled],||false).unwrap();
    let x=[0.3,-1.1,4.2];let h=1e-5;
    let mut points=vec![x];for c in 0..3 {
        let mut a=x;let mut b=x;a[c]+=h;b[c]-=h;points.push(a);points.push(b);
    }
    let rows=owner.receivers(&points,0.01).unwrap().prepare(0.5,medium(),options(),true,||false).unwrap();
    let work=rows.kernel_evaluations();let actual=rows.evaluate_velocity(&sources[0]).unwrap();
    let other=rows.evaluate_velocity(&sources[1]).unwrap();
    assert_eq!(work,rows.kernel_evaluations());
    for (a,b) in actual.scalar.pressure.iter().zip(&other.scalar.pressure) {
        assert!((*a*C64::new(-0.4,0.3)-*b).abs()<1e-10*(1.+a.abs()));
    }
    for c in 0..3 {
        let fd=(actual.scalar.pressure[1+2*c]-actual.scalar.pressure[2+2*c]).scale(0.5/h);
        let analytic=actual.particle_velocity_m_s[0][c]*C64::new(0.,0.5*343.*1.2);
        assert!((fd-analytic).abs()<3e-6*analytic.abs().max(1e-9));
    }
    let original=sources[0].pressure_at(x).unwrap();
    assert!((actual.scalar.pressure[0]-original).abs()<5e-5*original.abs());
    let scalar=owner.receivers(&[x],0.01).unwrap().prepare(0.5,medium(),options(),false,||false).unwrap();
    assert!((scalar.evaluate(&sources[0]).unwrap().pressure[0]-actual.scalar.pressure[0]).abs()<1e-7*original.abs());
    assert!(scalar.evaluate_velocity(&sources[0]).is_err());
}

#[test]
fn g4_geometry_source_frequency_and_refinement_refusals_never_poison_retry() {
    let mesh=surface();let owner=Prepared::new(&mesh,0.7,budget(),||false).unwrap();
    for x in [[0.1,0.1,0.1],[0.2,0.2,0.],[f64::NAN,0.,1.]] {
        assert!(owner.receivers(&[x],0.001).is_err());
    }
    let geometry=owner.receivers(&[[0.22,0.18,-0.08]],0.001).unwrap();
    assert!(geometry.prepare(0.8,medium(),options(),true,||false).is_err());
    assert!(geometry.prepare(0.7,medium(),Options {maximum_depth:0,..options()},true,||false).is_err());
    assert!(geometry.prepare(0.7,medium(),Options {maximum_kernel_evaluations:80,..options()},true,||false).is_err());
    let mut calls=0;assert!(matches!(geometry.prepare(0.7,medium(),options(),true,
        ||{calls+=1;calls==7}),Err(HelmholtzError::Cancelled)));
    let rows=geometry.prepare(0.7,medium(),options(),true,||false).unwrap();
    let mut source=manufactured(&owner);let before=source.coefficients.clone();
    let first=rows.evaluate_velocity(&source).unwrap();
    let other=Prepared::new(&mesh,0.7,budget(),||false).unwrap();
    assert!(rows.evaluate(&manufactured(&other)).is_err(),"identical dimensions are not source identity");
    source.k=0.6;assert!(rows.evaluate(&source).is_err());source.k=0.7;
    source.medium.density=1.;assert!(rows.evaluate(&source).is_err());source.medium=medium();
    source.coefficients[0]=C64::new(f64::NAN,0.);assert!(rows.evaluate(&source).is_err());
    source.coefficients.copy_from_slice(&before);
    let retry=geometry.prepare(0.7,medium(),options(),true,||false).unwrap();
    let after=retry.evaluate_velocity(&source).unwrap();
    assert_eq!(rows.kernel_evaluations(),retry.kernel_evaluations());
    assert_eq!(first.scalar.pressure,after.scalar.pressure);
    assert_eq!(first.particle_velocity_m_s,after.particle_velocity_m_s);
    assert_eq!(source.coefficients,before);
}
