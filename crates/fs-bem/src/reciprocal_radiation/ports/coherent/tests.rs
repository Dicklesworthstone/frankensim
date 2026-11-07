//! G0/G3: original coherent work, outgoing flux and adaptive receiver fields.
use super::*;
use crate::{
    helmholtz::Medium,
    near_field::{FirstOrder, Options},
    panel3d::SpherePanels,
    reciprocal_radiation::{Budget, Prepared},
};

fn budget() -> Budget {
    Budget { maximum_unknowns: 128, maximum_cache_bytes: 8 << 20,
        maximum_work: 1_000_000_000, angular_refinement: 0 }
}
fn limits() -> CoherentPortFieldLimits {
    CoherentPortFieldLimits { maximum_payload_bytes: 1 << 20, maximum_work: 10_000_000 }
}
fn medium() -> Medium { Medium { density: 1.0, sound_speed: 1.0 } }
fn rows(surface: &SpherePanels) -> Vec<Vec<f64>> {
    let first: Vec<_> = surface.normals().iter().map(|n| n[0] + 0.2).collect();
    vec![first.clone(), first.iter().map(|v| -0.7 * v).collect(),
        surface.normals().iter().map(|n| n[2]).collect(), vec![0.0; first.len()]]
}
fn near(a: C64, b: C64) {
    assert!((a-b).abs() < 3e-10 * (1.0 + b.abs()), "{a:?} != {b:?}");
}

#[test]
fn g3_coherent_densities_work_flux_and_receivers_match_one_direct_source_solve() {
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let prepared = Prepared::new(&surface, 0.8, budget(), || false).unwrap();
    let shape = rows(&surface);
    let refs: Vec<_> = shape.iter().map(Vec::as_slice).collect();
    let amplitudes = [C64::new(0.4,0.7), C64::new(-0.3,0.2),
        C64::new(0.8,-0.1), C64::new(2.0,1.0)];
    let coherent = {
        let source = prepared.solve_velocity_ports(0.6, medium(), &refs, || false).unwrap();
        source.coherent(&amplitudes, limits(), || false).unwrap()
    }; // Unit columns can be released; immutable geometry still owns the field.
    let velocity: Vec<_> = (0..surface.areas().len()).map(|p|
        (0..amplitudes.len()).fold(C64::ZERO, |sum,j| sum + amplitudes[j].scale(shape[j][p]))).collect();
    let direct = prepared.solve_batch(0.6, medium(), &[&velocity], || false).unwrap();
    for (&a,&b) in coherent.surface_pressure_pa().iter().zip(&direct[0].pressure) { near(a,b); }
    for (&a,&b) in coherent.surface_velocity_m_s().iter().zip(&velocity) { near(a,b); }
    near(coherent.surface_power_w(), direct[0].surface_power);
    let power = coherent.power(100_000_000, || false).unwrap();
    assert!((power.outgoing_power_w-direct[0].outgoing_power_w).abs()
        < 3e-10 * (1.0+direct[0].outgoing_power_w));
    assert!(power.relative_flux_defect < 1e-8);
    let geometry = prepared.receivers(&[[0.3,-0.5,2.0],[-2.0,0.2,0.3]], 0.01).unwrap();
    let observer = geometry.prepare(0.6, medium(), Options::default(), true, || false).unwrap();
    let work = observer.kernel_evaluations();
    let actual = observer.evaluate_coherent_velocity(&coherent).unwrap();
    let expected = observer.evaluate_velocity(&direct[0]).unwrap();
    for r in 0..2 {
        near(actual.scalar.pressure[r], expected.scalar.pressure[r]);
        for c in 0..3 { near(actual.particle_velocity_m_s[r][c], expected.particle_velocity_m_s[r][c]); }
    }
    assert_eq!(work, observer.kernel_evaluations());
    near(coherent.far_field([0.3,0.2,1.0]).unwrap(), direct[0].far_field([0.3,0.2,1.0]).unwrap());
}

#[test]
fn g3_coincident_directional_arrays_preserve_pressure_and_source_independence() {
    let surface = SpherePanels::icosphere(1.0, 0).unwrap();
    let prepared = Prepared::new(&surface, 0.8, budget(), || false).unwrap();
    let shape = rows(&surface);
    let refs: Vec<_> = shape.iter().map(Vec::as_slice).collect();
    let source = prepared.solve_velocity_ports(0.5, medium(), &refs, || false).unwrap();
    let original_impedance = source.impedance().to_vec();
    let field = source.coherent(&[C64::ONE,C64::ZERO,C64::new(0.2,-0.1),C64::ONE],
        limits(), || false).unwrap();
    let positions = vec![[0.3,-0.5,2.0]; crate::near_field::MAX_RECEIVERS];
    let geometry = prepared.receivers(&positions, 0.01).unwrap();
    let observer = geometry.prepare(0.5, medium(), Options::default(), true, || false).unwrap();
    let result = observer.evaluate_coherent_velocity(&field).unwrap();
    let single_geometry = prepared.receivers(&positions[..1], 0.01).unwrap();
    let single = single_geometry.prepare(0.5, medium(), Options::default(), false, || false).unwrap();
    let reference = single.evaluate_coherent(&field).unwrap().pressure[0];
    let front = FirstOrder::new(0.5,[0.0,0.0,1.0]).unwrap();
    let back = FirstOrder::new(0.5,[0.0,0.0,-1.0]).unwrap();
    for r in 0..positions.len() {
        near(result.scalar.pressure[r], reference);
        let a = front.observe(result.scalar.pressure[r], result.particle_velocity_m_s[r], medium()).unwrap();
        let b = back.observe(result.scalar.pressure[r], result.particle_velocity_m_s[r], medium()).unwrap();
        near(a+b, result.scalar.pressure[r]);
    }
    assert_eq!(source.impedance(), original_impedance.as_slice());
    assert!(prepared.receivers(&[],0.01).is_err());
    assert!(prepared.receivers(&vec![positions[0];65],0.01).is_err());
    let refused = geometry.prepare(0.5, medium(),
        Options { maximum_kernel_evaluations:80, ..Options::default() }, true, || false);
    assert!(refused.is_err());
}

#[test]
fn g4_budget_cancellation_and_cross_wired_receivers_refuse_without_mutating_sources() {
    let surface = SpherePanels::icosphere(1.0,0).unwrap();
    let prepared = Prepared::new(&surface,0.8,budget(),||false).unwrap();
    let shape = rows(&surface);
    let refs: Vec<_> = shape.iter().map(Vec::as_slice).collect();
    let source = prepared.solve_velocity_ports(0.5,medium(),&refs,||false).unwrap();
    let rates = [C64::ONE,C64::ZERO,C64::ZERO,C64::ONE];
    let before = source.impedance().to_vec();
    assert!(source.coherent(&rates[..3],limits(),||false).is_err());
    assert!(source.coherent(&rates,CoherentPortFieldLimits { maximum_payload_bytes:1,
        ..limits() },||false).is_err());
    let required = 32*prepared.unknowns() + 32*surface.areas().len() + 32*64 + 256;
    assert!(source.coherent(&rates,CoherentPortFieldLimits { maximum_payload_bytes:required-1,
        ..limits() },||false).is_err());
    assert!(source.coherent(&rates,CoherentPortFieldLimits { maximum_payload_bytes:required,
        ..limits() },||false).is_ok());
    let mut at = 0;
    assert!(matches!(source.coherent(&rates,limits(),||{at+=1;at==3}),Err(HelmholtzError::Cancelled)));
    let field = source.coherent(&rates,limits(),||false).unwrap();
    assert!(field.power(1,||false).is_err());
    assert!(matches!(field.power(100_000_000,||true),Err(HelmholtzError::Cancelled)));
    let other = Prepared::new(&surface,0.8,budget(),||false).unwrap();
    let geometry = other.receivers(&[[0.0,0.0,2.0]],0.01).unwrap();
    let observer = geometry.prepare(0.5,medium(),Options::default(),true,||false).unwrap();
    assert!(observer.evaluate_coherent(&field).is_err());
    assert!(observer.evaluate_coherent_velocity(&field).is_err());
    assert_eq!(source.impedance(),before.as_slice());
    let zero = source.coherent(&[C64::ZERO;4],limits(),||false).unwrap().power(100_000_000,||false).unwrap();
    assert_eq!(zero.surface_power_w,C64::ZERO);
    assert_eq!(zero.outgoing_power_w,0.0);
    assert_eq!(zero.relative_flux_defect,0.0);
}
