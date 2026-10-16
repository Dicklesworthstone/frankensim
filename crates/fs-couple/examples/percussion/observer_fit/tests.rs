use super::*;

fn grid(intervals: usize) -> Vec<f64> {
    (0..=4*intervals).map(|i| core::f64::consts::TAU*(40.0+1600.0*i as f64/(4*intervals) as f64)).collect()
}
fn source() -> DiscreteStateSpace {
    DiscreteStateSpace { n: 2, a: vec![0.95,0.0,0.0,0.7], b: vec![0.2,0.1],
        c: vec![0.3,-0.4], d: 0.1, e_leftover: 0.0, t_s: 1.0/48_000.0 }
}

#[test]
fn independent_audit_keeps_physical_phase_scale_and_direct_feedthrough() {
    let original=source(); let omega=grid(20);
    let values: Vec<_>=omega.iter().map(|&w| original.eval(w).unwrap().conj()).collect();
    let (f,r)=fit(&omega,&values,original.t_s,8).unwrap();
    assert_eq!(r.order,2); assert_eq!(r.attempts,1);
    assert!(r.selection_maximum<1e-5 && r.audit_maximum<1e-5 && r.audit_rms<1e-5);
    for w in [300.0,1900.0,7000.0] {
        assert!((f.eval(w).unwrap()-original.eval(w).unwrap()).abs()<1e-5);
    }
    let scaled: Vec<_>=values.iter().map(|v|v.scale(-3.0)).collect();
    let (g,s)=fit(&omega,&scaled,original.t_s,8).unwrap();
    assert_eq!(s.order,r.order);
    for &w in &omega { assert!((g.eval(w).unwrap()+f.eval(w).unwrap().scale(3.0)).abs()<1e-5); }
}

#[test]
fn unused_audit_points_cannot_change_order_poles_or_fit_normalization() {
    let original=source(); let omega=grid(20);
    let values: Vec<_>=omega.iter().map(|&w|original.eval(w).unwrap().conj()).collect();
    let mut poisoned=values.clone();
    for i in (1..poisoned.len()).step_by(2) { poisoned[i]=C64::new(100.0,-70.0); }
    let (a,ra,sa)=select(&omega,&values,original.t_s,8).unwrap();
    let (b,rb,sb)=select(&omega,&poisoned,original.t_s,8).unwrap();
    assert_eq!(ra.order,rb.order); assert_eq!(ra.attempts,rb.attempts); assert_eq!(sa,sb);
    assert_eq!(a.a,b.a); assert_eq!(a.b,b.b); assert_eq!(a.c,b.c); assert_eq!(a.d,b.d);
    assert!(fit(&omega,&poisoned,original.t_s,8).is_err());
}

#[test]
fn exact_zero_and_unseen_nonzero_sources_are_distinguished() {
    let omega=grid(20); let values=vec![C64::new(0.0,0.0);omega.len()];
    let (f,r)=fit(&omega,&values,1.0/48_000.0,8).unwrap();
    assert_eq!(f.n,0); assert_eq!(f.d,0.0); assert_eq!(r.order,0); assert_eq!(r.audit_maximum,0.0);
    let mut changed=values; changed[1]=C64::new(1e-200,0.0);
    assert!(fit(&omega,&changed,1.0/48_000.0,8).is_err());
}

#[test]
fn finite_data_grid_and_order_budgets_refuse_before_fitting() {
    let omega=grid(20); let values=vec![C64::new(1.0,0.0);omega.len()];
    for order in [0,1,3,12,34,usize::MAX] {
        assert!(fit(&omega,&values,1.0/48_000.0,order).is_err());
    }
    for dt in [0.0,-1.0,f64::INFINITY,1.0] { assert!(fit(&omega,&values,dt,8).is_err()); }
    let mut w=omega.clone(); w[3]=w[2]; assert!(fit(&w,&values,1.0/48_000.0,8).is_err());
    let mut v=values.clone(); v[1].re=f64::NAN; assert!(fit(&omega,&v,1.0/48_000.0,8).is_err());
    assert!(fit(&omega[..80],&values[..80],1.0/48_000.0,8).is_err());
    assert!(fit(&omega,&values[..80],1.0/48_000.0,8).is_err());
}
