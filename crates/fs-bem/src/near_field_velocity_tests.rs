use super::*;
use crate::helmholtz::{Formulation,solve_radiation};
fn triangles()->Vec<Triangle> {
    let p=[[0.,0.,-0.002],[0.1,0.,-0.002],[0.1,0.1,-0.002],[0.,0.1,-0.002],
        [0.,0.,0.002],[0.1,0.,0.002],[0.1,0.1,0.002],[0.,0.1,0.002]];
    [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
        [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]].map(|t|t.map(|i|p[i])).to_vec()
}
fn medium()->Medium {Medium {density:1.2,sound_speed:340.}}
fn solution(body:&SpherePanels,k:f64)->RadiationSolution {
    let velocity=body.normals().iter().map(|n|C64::new(n[2],0.)).collect::<Vec<_>>();
    solve_radiation(body,k,medium(),&velocity,Formulation::PlainCbie).unwrap()
}
#[test]
fn first_order_plane_wave_patterns_and_spherical_proximity_follow_euler() {
    let air=medium();let p=C64::new(1.3,-0.2);let velocity=[C64::ZERO,C64::ZERO,p.scale(1./(air.density*air.sound_speed))];
    let observe=|a,axis|FirstOrder::new(a,axis).unwrap().observe(p,velocity,air).unwrap();
    assert_eq!(observe(1.,[0.,0.,1.]),p);
    assert!((observe(0.5,[0.,0.,-1.])-p).abs()<1e-14);
    assert!(observe(0.5,[0.,0.,1.]).abs()<1e-14);
    assert!((observe(0.5,[1.,0.,0.])-p.scale(0.5)).abs()<1e-14);
    assert!((observe(0.,[0.,0.,1.])+p).abs()<1e-14);
    // Exact monopole Euler velocity: u=p/(rho*c)*(1+i/(kr)) along propagation.
    // The non-plane reactive term, not a distance gain, creates proximity.
    for kr in [0.02,0.2,2.,20.] {
        let v=[C64::ZERO,C64::ZERO,p*C64::new(1.,1./kr).scale(1./(air.density*air.sound_speed))];
        for a in [0.,0.5,1.] {
            let mic=FirstOrder::new(a,[0.,0.,-1.]).unwrap();
            let actual=mic.observe(p,v,air).unwrap();
            assert!((actual-p*C64::new(1.,(1.-a)/kr)).abs()<1e-12);
        }
    }
}
#[test]
fn analytic_velocity_matches_spatial_differences_of_the_same_bem_pressure() {
    let body=SpherePanels::from_triangles(triangles()).unwrap();let air=medium();let x=[0.031,0.043,0.027];
    let geometry=Geometry::new(&body,&[x],0.01).unwrap();let epsilon=1e-6;
    for k in [0.3,4.] {
        let source=solution(&body,k);let original=source.pressure.clone();
        let prepared=geometry.prepare_velocity(k,air,Options::default()).unwrap();
        let result=prepared.evaluate_velocity(&source).unwrap();
        let scalar=geometry.prepare(k,air,Options::default()).unwrap().evaluate(&source).unwrap();
        assert!((result.scalar.pressure[0]-scalar.pressure[0]).abs()<1e-7*scalar.pressure[0].abs());
        let measure=|point|Geometry::new(&body,&[point],0.01).unwrap().prepare(k,air,Options::default()).unwrap()
            .evaluate(&source).unwrap().pressure[0];
        for c in 0..3 {
            let mut a=x;let mut b=x;a[c]+=epsilon;b[c]-=epsilon;
            let fd=(measure(a)-measure(b))*C64::new(0.,-1./(2.*epsilon*k*air.sound_speed*air.density));
            let actual=result.particle_velocity_m_s[0][c];
            assert!((fd-actual).abs()<2e-6*actual.abs().max(1e-10),"axis {c}: {fd:?} vs {actual:?}");
            assert!(result.quadrature_error_estimate_m_s[0][c].is_finite());
        }
        let gain=C64::new(0.7,-0.3);let mut scaled=source.clone();
        for z in scaled.pressure.iter_mut().chain(&mut scaled.velocity) {*z=*z*gain;}
        let scaled=prepared.evaluate_velocity(&scaled).unwrap();
        for (a,b) in result.particle_velocity_m_s[0].iter().zip(scaled.particle_velocity_m_s[0]) {
            assert!((b-*a*gain).abs()<1e-12*(1.+a.abs()));
        }
        assert_eq!(source.pressure,original);
    }
}
#[test]
fn receiver_vector_and_oriented_response_are_rigid_frame_covariant() {
    let base=triangles();let body=SpherePanels::from_triangles(base.clone()).unwrap();let k=3.;let air=medium();
    let (s,c)=0.71_f64.sin_cos();let rot=|p:Point|[c*p[0]-s*p[2],p[1],s*p[0]+c*p[2]];
    let pose=|p:Point| {let v=rot(p);[v[0]+0.3,v[1]-0.2,v[2]+0.1]};
    let rotated=SpherePanels::from_triangles(base.iter().map(|t|t.map(pose)).collect()).unwrap();
    let source=solution(&body,k);let mut other=solution(&rotated,k);
    // Scalar pressure and normal-velocity traces are invariant under a proper
    // frame change. Only their known source-surface fingerprint changes.
    other.pressure.clone_from(&source.pressure);other.velocity.clone_from(&source.velocity);
    let point=[0.031,0.043,0.027];
    let a=Geometry::new(&body,&[point],0.01).unwrap().prepare_velocity(k,air,Options::default()).unwrap().evaluate_velocity(&source).unwrap();
    let b=Geometry::new(&rotated,&[pose(point)],0.01).unwrap().prepare_velocity(k,air,Options::default()).unwrap().evaluate_velocity(&other).unwrap();
    assert!((a.scalar.pressure[0]-b.scalar.pressure[0]).abs()<1e-7*a.scalar.pressure[0].abs());
    let expected=[a.particle_velocity_m_s[0][0].scale(c)-a.particle_velocity_m_s[0][2].scale(s),
        a.particle_velocity_m_s[0][1],a.particle_velocity_m_s[0][0].scale(s)+a.particle_velocity_m_s[0][2].scale(c)];
    for i in 0..3 {assert!((b.particle_velocity_m_s[0][i]-expected[i]).abs()<1e-7*expected[i].abs().max(1e-10));}
    let ma=FirstOrder::new(0.5,[0.,0.,-1.]).unwrap();let mb=FirstOrder::new(0.5,rot([0.,0.,-1.])).unwrap();
    let pa=ma.observe(a.scalar.pressure[0],a.particle_velocity_m_s[0],air).unwrap();
    let pb=mb.observe(b.scalar.pressure[0],b.particle_velocity_m_s[0],air).unwrap();
    assert!((pa-pb).abs()<1e-7*pa.abs());
}
#[test]
fn velocity_requires_explicit_preparation_and_preserves_identity_and_work_refusals() {
    let body=SpherePanels::from_triangles(triangles()).unwrap();let k=3.;let air=medium();
    let geometry=Geometry::new(&body,&[[0.04,0.04,0.022]],0.01).unwrap();let source=solution(&body,k);
    let scalar=geometry.prepare(k,air,Options::default()).unwrap();assert!(scalar.evaluate_velocity(&source).is_err());
    let prepared=geometry.prepare_velocity(k,air,Options::default()).unwrap();
    let mut bad_source=source.clone();bad_source.k*=2.;assert!(prepared.evaluate_velocity(&bad_source).is_err());
    bad_source=source.clone();bad_source.surface_fingerprint^=1;assert!(prepared.evaluate_velocity(&bad_source).is_err());
    bad_source=source;bad_source.velocity[0]=C64::new(f64::NAN,0.);assert!(prepared.evaluate_velocity(&bad_source).is_err());
    assert!(geometry.prepare_velocity(k,air,Options {maximum_kernel_evaluations:80,..Options::default()}).is_err());
    for a in [-1.,1.01,f64::NAN] {assert!(FirstOrder::new(a,[0.,0.,1.]).is_err());}
    for axis in [[0.;3],[0.,0.,2.],[f64::NAN,0.,1.]] {assert!(FirstOrder::new(0.5,axis).is_err());}
}
