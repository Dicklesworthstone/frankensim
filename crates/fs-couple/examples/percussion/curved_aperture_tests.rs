use super::*;
use fs_plate::shell::{head::TensionedDiskSpec,profile::ProfileBudget};

fn fixture() -> Boundary {
    let film=|| TensionedDisk::new(TensionedDiskSpec {radius_m:0.1,thickness_m:0.0002,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:1000.0,
        radial_intervals:2,azimuths:8},ProfileBudget {max_nodes:100,max_triangles:200,
        max_feature_evaluations:0}).unwrap();
    let films=[film(),film()];
    let modes: Vec<Vec<ModePair>>=films.iter().map(|f| {
        let mut phi=vec![0.0;f.model.free];
        for (i,&(x,y)) in f.mesh.nodes.iter().enumerate() {
            if let Some(k)=f.model.dof_map[3*i] {phi[k]=1.0-(x*x+y*y)/0.01;}
        }
        // Kinematic projection fixture, not an eigenfrequency claim.
        vec![ModePair {lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)}]
    }).collect();
    Boundary::drum(&films,&modes,0.12,0.11).unwrap()
}
fn opening(angle:f64) -> SidewallAperture {
    SidewallAperture {radius_m:0.005,azimuth_rad:angle,axial_position_m:0.067,
        radial_rings:8,angular_points:32,maximum_terms:2_000_000}
}
fn flux(b:&Boundary, mode:usize) -> f64 {
    b.triangles.iter().zip(&b.weights[mode]).map(|(t,w)|area(*t)*w).sum()
}
fn port() -> NeckRadiationPort {
    NeckRadiationPort {coordinate:7,area_m2:std::f64::consts::PI*0.005_f64.powi(2),
        effective_length_m:0.008,volume_weight_m2_per_sqrt_kg:0.003}
}
fn prepare(angle:f64) -> Boundary {
    fixture().with_unrolled_sidewall_aperture(0.11,0.12,opening(angle),port(),2048,
        &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn finite_mouth_conserves_neck_flow_and_both_original_head_projections() {
    let original=fixture(); let b=prepare(0.2);
    assert!(b.triangles.len()>original.triangles.len() && b.triangles.len()<=2048);
    assert_eq!(b.state_modes(), &[1,2,7]);
    for head in 0..2 {assert!((flux(&b,head)-flux(&original,head)).abs()<1e-14);}
    assert!((flux(&b,2)-0.003).abs()<1e-16);
    let mut sources=0;
    for (t,&w) in b.triangles.iter().zip(&b.weights[2]) {
        assert!(area(*t)>0.0);
        if w>0.0 {
            sources+=1;
            // Flux belongs to the outer wall, never a head or the air-pressure state.
            assert!(t.iter().all(|p| p[2].abs()<0.06));
            assert!(t.iter().all(|p| p[0].hypot(p[1])>0.1));
        }
    }
    assert!(sources>8,"the mouth must be spatially resolved, not snapped to one panel");
    // A negative actual slug velocity reverses the full physical flux.
    assert!((flux(&b,2)*-2.0+0.006).abs()<2e-16);
}

#[test]
fn seam_crossing_and_periodic_angles_keep_a_local_conservative_source() {
    let angle=std::f64::consts::PI-0.006;let b=prepare(angle);
    let other=prepare(angle-std::f64::consts::TAU);
    assert_eq!(b.triangles.len(),other.triangles.len());
    assert!((flux(&b,2)-flux(&other,2)).abs()<1e-16);
    let mut center=[0.0;3];
    for (t,&w) in b.triangles.iter().zip(&b.weights[2]) {
        for axis in 0..3 {center[axis]+=area(*t)*w*t.iter().map(|p|p[axis]/3.0).sum::<f64>()/0.003;}
    }
    let expected=[0.11*angle.cos(),0.11*angle.sin(),0.06-0.067];
    for axis in 0..3 {assert!((center[axis]-expected[axis]).abs()<0.004);}
    // The first head's outward sign remains negative and the second positive.
    assert!(flux(&b,0)<0.0 && flux(&b,1)>0.0);
}

#[test]
fn malformed_missing_and_over_budget_apertures_refuse_without_fallback() {
    let gate=CancelGate::new_clock_free();
    for a in [SidewallAperture {radius_m:0.0,..opening(0.0)},
        SidewallAperture {radius_m:0.02,..opening(0.0)},
        SidewallAperture {axial_position_m:0.001,..opening(0.0)},
        SidewallAperture {azimuth_rad:f64::NAN,..opening(0.0)},
        SidewallAperture {maximum_terms:1,..opening(0.0)}] {
        assert!(fixture().with_unrolled_sidewall_aperture(0.11,0.12,a,port(),2048,&gate).is_err());
    }
    for bad in [NeckRadiationPort {coordinate:1,..port()},
        NeckRadiationPort {volume_weight_m2_per_sqrt_kg:0.0,..port()},
        NeckRadiationPort {volume_weight_m2_per_sqrt_kg:f64::NAN,..port()},
        NeckRadiationPort {coordinate:usize::MAX,..port()},
        NeckRadiationPort {effective_length_m:0.02,..port()}] {
        assert!(fixture().with_unrolled_sidewall_aperture(0.11,0.12,opening(0.0),bad,2048,&gate).is_err());
    }
    let original=fixture();let count=original.triangles.len();
    assert!(original.with_unrolled_sidewall_aperture(0.11,0.12,opening(0.0),port(),count,&gate).is_err());
    // A wrong cylinder radius must not move/snap the source onto some other surface.
    assert!(fixture().with_unrolled_sidewall_aperture(0.115,0.12,opening(0.0),port(),2048,&gate).is_err());
    gate.request();assert!(fixture().with_unrolled_sidewall_aperture(0.11,0.12,opening(0.0),port(),2048,&gate).is_err());
}

#[test]
fn aperture_flow_reaches_the_real_exterior_neumann_pressure_operator() {
    let b=prepare(0.2);let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let omega=std::f64::consts::TAU*90.0;let medium=Medium::air();
    let positive=acceleration_fields(&b.weights[2..],omega).remove(0);
    let negative:Vec<_>=positive.iter().map(|p|p.scale(-2.0)).collect();
    let solutions=solve_radiation_batch(&surface,omega/medium.sound_speed,medium,
        &[&positive,&negative],Formulation::PlainCbie).unwrap();
    let receiver=[[0.3,0.0,0.1]];
    let p=exterior_pressure_at_points(&surface,&solutions[0],medium,&receiver).unwrap()[0];
    let reversed=exterior_pressure_at_points(&surface,&solutions[1],medium,&receiver).unwrap()[0];
    assert!(p.abs()>1e-10 && p.abs().is_finite());
    assert!((reversed+p.scale(2.0)).abs()<1e-8*p.abs());
    assert!(solutions[0].radiated_power_roundoff_interval.1>=0.0);
}
