use super::*;
use crate::{panel3d::SpherePanels, reciprocal_radiation::Budget};
fn budget() -> Budget { Budget {maximum_unknowns:128,maximum_cache_bytes:8*1024*1024,
    maximum_work:1_000_000_000,angular_refinement:0} }
fn medium() -> Medium { Medium {density:1.0,sound_speed:1.0} }
fn rows(surface: &SpherePanels) -> Vec<Vec<f64>> {
    let a: Vec<_> = surface.normals().iter().map(|n| n[0]+0.2).collect();
    vec![a.clone(),a.iter().map(|v|-0.7*v).collect(),
        surface.normals().iter().map(|n|n[2]).collect(),vec![0.0;a.len()]]
}
#[test]
fn g0_ports_preserve_signed_dependent_and_zero_columns_with_original_solution() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let prepared=Prepared::new(&surface,0.8,budget(),||false).unwrap();
    let rows=rows(&surface);let refs:Vec<_>=rows.iter().map(Vec::as_slice).collect();
    let result=prepared.solve_velocity_ports(0.5,medium(),&refs,||false).unwrap();
    assert_eq!(result.ports(),4);let z=result.impedance();
    for i in 0..4 {for j in 0..4 {
        let manual: C64=rows[i].iter().zip(&result.solutions()[j].pressure).zip(surface.areas())
            .fold(C64::ZERO,|s,((b,p),a)|s+p.scale(b*a));
        assert_eq!(z[4*i+j],manual);
        assert!((z[4*i+j]-z[4*j+i]).abs()<1e-10*(1.0+z[4*i+j].abs()));
    }}
    assert!(z[1].re<0.0);
    for i in 0..4 {assert_eq!(z[4*i+3],C64::ZERO);assert_eq!(z[12+i],C64::ZERO);}
    let receivers=prepared.receivers(&[[0.0,0.0,2.0]],1e-8).unwrap();
    let observation=receivers.prepare(0.5,medium(),crate::near_field::Options::default(),false,||false).unwrap();
    let first=observation.evaluate(&result.solutions()[0]).unwrap().pressure[0];
    let second=observation.evaluate(&result.solutions()[1]).unwrap().pressure[0];
    assert!((second-first.scale(-0.7)).abs()<1e-10*(1.0+first.abs()));
}
#[test]
fn g3_complex_port_mixtures_reproduce_actual_surface_work() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let prepared=Prepared::new(&surface,0.8,budget(),||false).unwrap();
    let rows=rows(&surface);let refs:Vec<_>=rows.iter().map(Vec::as_slice).collect();
    let result=prepared.solve_velocity_ports(0.6,medium(),&refs,||false).unwrap();
    let amplitudes=[C64::new(0.4,0.7),C64::new(-0.3,0.2),C64::new(0.8,-0.1),C64::new(2.0,1.0)];
    let velocity:Vec<_>=(0..surface.areas().len()).map(|p|(0..4).fold(C64::ZERO,
        |v,j|v+amplitudes[j].scale(rows[j][p]))).collect();
    let direct=prepared.solve_batch(0.6,medium(),&[&velocity],||false).unwrap();
    let mut work=C64::ZERO;
    for i in 0..4 {for j in 0..4 {
        work=work+amplitudes[i].conj()*result.impedance()[4*i+j]*amplitudes[j];
    }}work=work.scale(0.5);
    assert!((work-direct[0].surface_power).abs()<1e-10*(1.0+work.abs()));
    assert!(work.re>=0.0);
    let receiver=[0.3,-0.5,2.0];
    let pressure=result.solutions().iter().zip(amplitudes).fold(C64::ZERO,
        |p,(s,a)|p+s.pressure_at(receiver).unwrap()*a);
    assert!((pressure-direct[0].pressure_at(receiver).unwrap()).abs()<1e-10*(1.0+pressure.abs()));
}
#[test]
fn g4_bad_rows_and_cancelled_projection_do_not_poison_retry() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let prepared=Prepared::new(&surface,0.8,budget(),||false).unwrap();
    let rows=rows(&surface);let refs:Vec<_>=rows.iter().map(Vec::as_slice).collect();
    assert!(prepared.solve_velocity_ports(0.5,medium(),&[],||false).is_err());
    assert!(prepared.solve_velocity_ports(0.5,medium(),&[&[0.0]],||false).is_err());
    let mut invalid=rows[0].clone();invalid[1]=f64::NAN;
    assert!(prepared.solve_velocity_ports(0.5,medium(),&[&invalid],||false).is_err());
    let mut polls=0;
    let first=prepared.solve_velocity_ports(0.5,medium(),&refs,||{polls+=1;false}).unwrap();
    let mut at=0;
    assert!(matches!(prepared.solve_velocity_ports(0.5,medium(),&refs,||{
        at+=1;at==polls-1 // after the source solve, during final work projection
    }),Err(HelmholtzError::Cancelled)));
    let retry=prepared.solve_velocity_ports(0.5,medium(),&refs,||false).unwrap();
    assert_eq!(first.impedance(),retry.impedance());
    assert_eq!(first.solutions()[0].pressure,retry.solutions()[0].pressure);
}
#[test]
fn g4_projection_work_is_included_in_the_admitted_frequency_budget() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let prepared=Prepared::new(&surface,0.8,budget(),||false).unwrap();
    let solve=prepared.frequency_work(4).unwrap();
    assert!(prepared.velocity_port_work(4).unwrap()>solve);
    let constrained=Prepared::new(&surface,0.8,Budget{maximum_work:solve,..budget()},||false).unwrap();
    let rows=rows(&surface);let refs:Vec<_>=rows.iter().map(Vec::as_slice).collect();
    assert!(constrained.solve_velocity_ports(0.5,medium(),&refs,||false).is_err());
}
