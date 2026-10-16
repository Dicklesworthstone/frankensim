use fs_material::sheet::{SheetLaw,QuadraticMembrane};
use fs_material::sheet::tabulated::{StrainCurve,TabulatedMembrane};
fn curve(k:f64)->StrainCurve {
    StrainCurve::new(&[(-0.2,-0.2*k),(0.0,0.0),(0.2,0.2*k)]).unwrap()
}
fn close(a:f64,b:f64) {assert!((a-b).abs()<1e-9+1e-9*a.abs().max(b.abs()),"{a} != {b}");}
#[test]
fn linear_tables_reproduce_the_original_quadratic_law() {
    let q=QuadraticMembrane::new(8000.0,2000.0,0.2,500.0).unwrap();
    let d=q.matrix().coefficients();
    let table=TabulatedMembrane::new([curve(d[0]),curve(d[4]),curve(d[8])],d[1]).unwrap();
    for strain in [[0.0;3],[-0.03,0.01,0.04],[0.1,-0.08,-0.02]] {
        let a=q.response(strain).unwrap();let b=table.response(strain).unwrap();
        close(a.energy_j_m2,b.energy_j_m2);
        for (a,b) in a.resultants_n_m.iter().zip(b.resultants_n_m) {close(*a,b);}
        for (a,b) in a.tangent_n_m.iter().zip(b.tangent_n_m) {close(*a,b);}
    }
}
#[test]
fn nonlinear_energy_stress_and_tangent_share_one_curve() {
    let warp=StrainCurve::new(&[(-0.2,-80.0),(0.0,0.0),(0.02,20.0),(0.1,420.0)]).unwrap();
    let law=TabulatedMembrane::new([warp,curve(1200.0),curve(600.0)],50.0).unwrap();
    for x in [-0.08,0.01,0.04] {
        let e=[x,0.013,-0.007];let r=law.response(e).unwrap();let h=1e-7;
        for j in 0..3 {
            let mut lo=e;let mut hi=e;lo[j]-=h;hi[j]+=h;
            let a=law.response(lo).unwrap();let b=law.response(hi).unwrap();
            close((b.energy_j_m2-a.energy_j_m2)/(2.0*h),r.resultants_n_m[j]);
            for i in 0..3 {close((b.resultants_n_m[i]-a.resultants_n_m[i])/(2.0*h),r.tangent_n_m[3*i+j]);}
        }
    }
    assert_eq!(law.response([0.01,0.0,0.0]).unwrap().tangent_n_m[0],1000.0);
    assert_eq!(law.response([0.04,0.0,0.0]).unwrap().tangent_n_m[0],5000.0);
}
#[test]
fn origin_and_knot_semantics_preserve_work_without_cancellation() {
    let c=StrainCurve::new(&[(-1.0,-1000.0),(0.0,0.0),(0.01,10.0),(0.02,40.0)]).unwrap();
    assert_eq!(c.evaluate(0.0).unwrap(),[0.0,0.0,1000.0]);
    assert_eq!(c.evaluate(0.01).unwrap(),[0.05,10.0,3000.0]);
    assert_eq!(c.evaluate(0.02).unwrap(),[0.3,40.0,3000.0]);
    for e in [-1e-12,1e-12] {
        let r=c.evaluate(e).unwrap();assert!(r[0]>0.0);
        assert!((r[0]/(500.0*e*e)-1.0).abs()<1e-14);
    }
}
#[test]
fn unsupported_data_domain_and_unstable_coupling_refuse() {
    for points in [vec![(0.0,0.0),(1.0,1.0)],vec![(-1.0,-1.0),(0.0,0.0),(0.0,1.0)],
        vec![(-1.0,-1.0),(0.0,0.0),(1.0,0.0)],vec![(-1.0,-1.0),(0.0,1.0),(1.0,2.0)]] {
        assert!(StrainCurve::new(&points).is_err());
    }
    let law=TabulatedMembrane::new([curve(10.0),curve(10.0),curve(10.0)],0.0).unwrap();
    for e in [[0.21,0.0,0.0],[0.0,-0.21,0.0],[0.0,0.0,f64::NAN]] {assert!(law.response(e).is_err());}
    assert!(TabulatedMembrane::new([curve(10.0),curve(10.0),curve(10.0)],10.0).is_err());
}

#[test]
fn one_sided_newton_tangents_do_not_become_unsupported_vibration_frequencies() {
    let curve=StrainCurve::new(&[(-0.1,-100.0),(0.0,0.0),(0.02,20.0),(0.1,180.0)]).unwrap();
    assert_eq!(curve.evaluate(0.02).unwrap()[2],2000.0);
    assert!(curve.admit_linearization(0.02).is_err());
    assert!(curve.admit_linearization(0.02+f64::EPSILON).is_err());
    assert!(curve.admit_linearization(-0.1).is_err());
    assert!(curve.admit_linearization(0.1).is_err());
    assert!(curve.admit_linearization(0.0).is_ok()); // equal slopes across origin
    assert!(curve.admit_linearization(0.021).is_ok());
}
