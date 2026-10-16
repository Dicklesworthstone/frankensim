//! Geometry/basis regressions; analytic roots are independent reference values.
use fs_couple::render::plate::impact::{ImpactError,cavity::cylinder::{CylinderSpec,CylindricalCavity,SidewallAperture}};
use fs_couple::vibroacoustic::AcousticMedium;
use fs_exec::CancelGate;
fn spec(intervals:usize)->CylinderSpec {CylinderSpec {radius_m:0.1703,depth_m:0.1651,
    radial_intervals:intervals,maximum_azimuthal_order:2,maximum_axial_order:1,
    maximum_frequency_hz:1400.0,maximum_modes:32,eigen_residual_tolerance:1e-7}}
fn build(s:CylinderSpec)->CylindricalCavity {
    CylindricalCavity::new(s,AcousticMedium {rho0:1.2,c0:343.0},&CancelGate::new_clock_free()).unwrap()
}
#[test]
fn radial_frequencies_converge_to_neumann_bessel_roots_not_dirichlet_roots() {
    // First positive J_m' zeros, NIST DLMF 10.21, not fitted drum frequencies.
    let roots=[3.8317059702075125,1.8411837813406595,3.0542369282271404];
    let mut previous=[f64::INFINITY;3];
    for intervals in [8,16,32] {
        let air=build(spec(intervals));
        for order in 0..3 {
            let mode=air.modes().iter().find(|m|m.azimuthal_order==order && m.axial_order==0
                && !m.sine && m.radial_index==usize::from(order==0)).unwrap();
            let root=mode.omega_rad_s*spec(intervals).radius_m/343.0;
            let error=(root/roots[order]-1.0).abs();
            assert!(error<0.007 && error<0.35*previous[order]);
            if intervals==32 {assert!(error<0.0004);}
            previous[order]=error;
        }
    }
}
#[test]
fn uniform_compliance_axial_parity_and_angular_pairs_keep_physical_normalization() {
    let air=build(spec(32));let s=air.spec();
    let uniform=&air.modes()[0];
    let volume=core::f64::consts::PI*s.radius_m*s.radius_m*s.depth_m;
    assert_eq!(uniform.omega_rad_s,0.0);assert!((uniform.norm_m3-volume).abs()<1e-16);
    let top=air.values_at([0.08,0.04,0.0]).unwrap();
    let bottom=air.values_at([0.08,0.04,s.depth_m]).unwrap();
    assert_eq!(top[0],1.0);assert_eq!(bottom[0],1.0);
    let axial=air.modes().iter().position(|m|m.azimuthal_order==0 && m.radial_index==0 && m.axial_order==1).unwrap();
    assert!((air.modes()[axial].omega_rad_s-core::f64::consts::PI*343.0/s.depth_m).abs()<1e-10);
    assert!((air.modes()[axial].norm_m3-volume/2.0).abs()<1e-16);
    assert!((top[axial]+bottom[axial]).abs()<1e-12);
    for order in [1,2] {
        let i=air.modes().iter().position(|m|m.azimuthal_order==order && m.axial_order==0 && m.radial_index==0 && !m.sine).unwrap();
        assert!(air.modes()[i+1].sine);
        assert_eq!(air.modes()[i].omega_rad_s.to_bits(),air.modes()[i+1].omega_rad_s.to_bits());
        assert_eq!(air.modes()[i].norm_m3.to_bits(),air.modes()[i+1].norm_m3.to_bits());
        let rotated=air.values_at([-0.04,0.08,0.0]).unwrap();
        assert!((top[i]*top[i]+top[i+1]*top[i+1]-rotated[i]*rotated[i]-rotated[i+1]*rotated[i+1]).abs()<1e-12);
        assert_eq!(air.values_at([0.0,0.0,0.0]).unwrap()[i],0.0);
    }
}
#[test]
fn scaling_geometry_and_sound_speed_changes_derived_frequencies_and_norms() {
    let a=build(spec(16));let mut large=spec(16);large.radius_m*=2.0;large.depth_m*=2.0;large.maximum_frequency_hz/=2.0;
    let b=build(large);assert_eq!(a.modes().len(),b.modes().len());
    let av=a.values_at([0.06,0.03,0.08]).unwrap();let bv=b.values_at([0.12,0.06,0.16]).unwrap();
    for ((x,y),(u,v)) in a.modes().iter().zip(b.modes()).zip(av.iter().zip(&bv)) {
        assert!((x.omega_rad_s-2.0*y.omega_rad_s).abs()<1e-9);
        assert!((8.0*x.norm_m3-y.norm_m3).abs()<1e-14);assert!((u-v).abs()<1e-12);
    }
    let faster=CylindricalCavity::new(CylinderSpec {maximum_frequency_hz:2800.0,..spec(16)},
        AcousticMedium {rho0:1.2,c0:686.0},&CancelGate::new_clock_free()).unwrap();
    for (x,y) in a.modes().iter().zip(faster.modes()) {assert_eq!((2.0*x.omega_rad_s).to_bits(),y.omega_rad_s.to_bits());}
}
#[test]
fn complete_pair_capacity_bad_geometry_points_and_cancellation_refuse_explicitly() {
    let mut s=spec(16);s.maximum_modes=2;
    assert!(build_result(s).is_err()); // Cannot fit constant + both m=1 members.
    for radius in [0.0,-1.0,f64::NAN] {assert!(build_result(CylinderSpec {radius_m:radius,..spec(16)}).is_err());}
    let air=build(spec(16));
    assert!(air.values_at([air.spec().radius_m*1.001,0.0,0.0]).is_err());
    assert!(air.values_at([0.0,0.0,-0.1]).is_err());
    assert!(air.sample(&[[0.0,0.0,0.0]],0).is_err());
    let sample=air.sample(&[[0.0,0.0,0.0],[0.08,0.01,0.08]],100).unwrap();
    assert_eq!(sample.interface[0],[1.0,1.0]);assert_eq!(sample.loss_factor,0.0);
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(matches!(CylindricalCavity::new(spec(16),AcousticMedium {rho0:1.2,c0:343.0},&gate),Err(ImpactError::Cancelled)));
}
fn build_result(s:CylinderSpec)->Result<CylindricalCavity,ImpactError> {
    CylindricalCavity::new(s,AcousticMedium {rho0:1.2,c0:343.0},&CancelGate::new_clock_free())
}

fn opening()->SidewallAperture {SidewallAperture {radius_m:0.016,azimuth_rad:0.4,
    axial_position_m:0.06,radial_rings:8,angular_points:32,maximum_terms:100000}}

#[test]
fn g1_finite_sidewall_area_converges_to_the_independent_fourier_disk_integral() {
    let air=build(spec(16));let gate=CancelGate::new_clock_free();let s=air.spec();let aperture=opening();
    let centre=air.values_at([s.radius_m*aperture.azimuth_rad.cos(),
        s.radius_m*aperture.azimuth_rad.sin(),aperture.axial_position_m]).unwrap();
    // The disk characteristic function is 2 J1(x)/x. Independent convergent
    // reference series, NOT a second implementation of the production quadrature.
    let exact:Vec<_>=air.modes().iter().zip(&centre).map(|(mode,value)| {
        let x=aperture.radius_m*(mode.azimuthal_order as f64/s.radius_m)
            .hypot(core::f64::consts::PI*mode.axial_order as f64/s.depth_m);
        let mut term=1.0;let mut factor=1.0;
        for k in 1..=18 {term*= -x*x/(4.0*f64::from(k)*f64::from(k+1));factor+=term;}
        value*factor
    }).collect();
    let mut previous=f64::INFINITY;
    for radial_rings in [2,4,8,16] {
        let mean=air.sidewall_averages(SidewallAperture {radial_rings,..aperture},&gate).unwrap();
        assert_eq!(mean[0],1.0,"constant pressure gives exactly the physical opening area");
        let error=mean.iter().zip(&exact).map(|(a,b)|(a-b).abs()).fold(0.0_f64,f64::max);
        assert!(error<0.35*previous);previous=error;
    }
    assert!(previous<1e-7);
    assert!(exact.iter().zip(&centre).any(|(a,b)|(a-b).abs()>1e-3),"a finite opening cannot be replaced by a point");
    let original=air.sidewall_averages(aperture,&gate).unwrap();
    let opposite=air.sidewall_averages(SidewallAperture {
        azimuth_rad:aperture.azimuth_rad+core::f64::consts::PI,
        axial_position_m:s.depth_m-aperture.axial_position_m,..aperture},&gate).unwrap();
    for ((a,b),mode) in original.iter().zip(opposite).zip(air.modes()) {
        let parity=if (mode.azimuthal_order+mode.axial_order)%2==0 {1.0}else{-1.0};
        assert!((a-parity*b).abs()<1e-12);
    }
}

#[test]
fn g0_g4_sidewall_average_refuses_outside_geometry_work_and_cancellation() {
    let air=build(spec(8));let gate=CancelGate::new_clock_free();let good=opening();
    for bad in [SidewallAperture {radius_m:0.0,..good},
        SidewallAperture {radius_m:0.02,..good},SidewallAperture {axial_position_m:0.0,..good},
        SidewallAperture {azimuth_rad:f64::NAN,..good},SidewallAperture {radial_rings:0,..good},
        SidewallAperture {angular_points:7,..good},SidewallAperture {maximum_terms:1,..good}] {
        assert!(air.sidewall_averages(bad,&gate).is_err());
    }
    let expected=air.sidewall_averages(good,&gate).unwrap();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(air.sidewall_averages(good,&cancelled),Err(ImpactError::Cancelled)));
    assert_eq!(expected,air.sidewall_averages(good,&gate).unwrap());
}
