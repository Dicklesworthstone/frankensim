use super::*;
use crate::render::plate::impact::{BodyPotential,ImpactBody,ImpactConfig};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::plate::impact::cavity::CavityCoupling;
fn budget(n:usize)->Budget {Budget {radial_intervals:n,maximum_azimuthal_order:2,
    maximum_frequency_hz:2.0,maximum_modes:128,maximum_nodes:512,
    maximum_matrix_entries:1_000_000,relative_tolerance:1e-7}}
fn medium()->AcousticMedium {AcousticMedium {rho0:1.2,c0:1.0}}
fn cylinder(n:usize)->Vec<Station> {(0..=n).map(|i|Station {z_m:i as f64/n as f64,radius_m:1.0}).collect()}
fn cavity(p:Vec<Station>,b:Budget)->MeridianCavity {MeridianCavity::new(p,medium(),b,&CancelGate::new_clock_free()).unwrap()}

#[test]
fn g1_meridional_cylinder_refines_to_independent_neumann_roots_and_keeps_volume() {
    // Exact axial pi and first positive J_1'/J_2' zeros; not another FEM pencil.
    let roots=[PI,1.8411837813406595,3.0542369282271404];
    let mut previous=[f64::INFINITY;3];
    for n in [4,8,12] {
        let air=cavity(cylinder(n),budget(n));
        assert!((air.volume_m3()-PI).abs()<1e-13);
        assert_eq!(air.modes()[0].omega_rad_s,0.0);
        assert_eq!(air.modes()[0].norm_m3,air.volume_m3());
        for order in 0..3 {
            let mode=air.modes().iter().find(|m|m.azimuthal_order==order
                && m.family_index==usize::from(order==0) && !m.sine).unwrap();
            let error=(mode.omega_rad_s/roots[order]-1.0).abs();
            assert!(error<0.025 && error<0.55*previous[order],"order={order}, n={n}, error={error}");
            previous[order]=error;
        }
    }
}
#[test]
fn g3_profile_scaling_and_rotation_preserve_the_same_shapes_and_complete_pairs() {
    let profile=vec![Station {z_m:-1.,radius_m:0.},Station {z_m:-0.5,radius_m:0.7},Station {z_m:0.,radius_m:1.}];
    let a=cavity(profile.clone(),budget(4));
    let large=profile.iter().map(|p|Station {z_m:2.*p.z_m+3.,radius_m:2.*p.radius_m}).collect();
    let b=cavity(large,Budget {maximum_frequency_hz:1.,..budget(4)});
    assert_eq!(a.modes.len(),b.modes.len());
    let av=a.values_at([0.2,0.1,-0.25]).unwrap();let bv=b.values_at([0.4,0.2,2.5]).unwrap();
    let rotated=a.values_at([-0.1,0.2,-0.25]).unwrap();
    for (i,(x,y)) in a.modes.iter().zip(&b.modes).enumerate() {
        assert!((x.omega_rad_s-2.*y.omega_rad_s).abs()<1e-8);
        assert!((8.*x.norm_m3-y.norm_m3).abs()<1e-10);
        assert!((av[i]-bv[i]).abs()<1e-9);
        if x.azimuthal_order>0 && !x.sine {
            let next=&a.modes[i+1];assert!(next.sine);
            assert_eq!(x.omega_rad_s.to_bits(),next.omega_rad_s.to_bits());
            assert_eq!(x.norm_m3.to_bits(),next.norm_m3.to_bits());
            assert!((av[i]*av[i]+av[i+1]*av[i+1]-rotated[i]*rotated[i]-rotated[i+1]*rotated[i+1]).abs()<1e-10);
        }
    }
    assert_eq!(a.values_at([0.,0.,-1.]).unwrap()[0],1.);
    assert!(a.values_at([0.9,0.,-0.8]).is_err(),"inside bounding cylinder but outside the bowl");
    assert!(a.values_at([0.,0.,0.01]).is_err());
}
#[test]
fn g3_same_height_and_lip_do_not_replace_bowl_volume_or_spectrum_with_a_cylinder() {
    let bowl=cavity(vec![Station {z_m:0.,radius_m:0.},Station {z_m:0.5,radius_m:0.7},
        Station {z_m:1.,radius_m:1.}],budget(8));
    let pipe=cavity(vec![Station {z_m:0.,radius_m:1.},Station {z_m:0.5,radius_m:1.},
        Station {z_m:1.,radius_m:1.}],budget(8));
    assert!(bowl.volume_m3()<0.6*pipe.volume_m3());
    let first=|air:&MeridianCavity|air.modes.iter().find(|m|m.azimuthal_order==1).unwrap().omega_rad_s;
    assert!((first(&bowl)/first(&pipe)-1.).abs()>0.1);
}
#[test]
fn g4_invalid_profiles_incomplete_pairs_and_cancelled_sampling_return_no_model() {
    let g=CancelGate::new_clock_free();
    for p in [vec![],vec![Station {z_m:0.,radius_m:1.}],
        vec![Station {z_m:0.,radius_m:1.};2],
        vec![Station {z_m:0.,radius_m:1.},Station {z_m:1.,radius_m:-1.}],
        vec![Station {z_m:0.,radius_m:1.},Station {z_m:0.5,radius_m:0.},Station {z_m:1.,radius_m:1.}]] {
        assert!(MeridianCavity::new(p,medium(),budget(4),&g).is_err());
    }
    for b in [Budget {maximum_nodes:4,..budget(4)},Budget {maximum_matrix_entries:1,..budget(4)},
        Budget {maximum_frequency_hz:0.4,maximum_modes:2,..budget(4)}] {
        assert!(MeridianCavity::new(cylinder(4),medium(),b,&g).is_err());
    }
    let air=cavity(cylinder(4),budget(4));let points=[[0.1,0.2,0.3]];
    let expected=air.sample(&points,128,&g).unwrap();assert_eq!(expected.interface[0],[1.]);
    assert!(air.sample(&points,0,&g).is_err());
    g.request();assert!(matches!(air.sample(&points,128,&g),Err(ImpactError::Cancelled)));
    assert!(matches!(MeridianCavity::new(cylinder(4),medium(),budget(4),&g),Err(ImpactError::Cancelled)));
    let actual=air.sample(&points,128,&CancelGate::new_clock_free()).unwrap();
    assert_eq!(actual.interface,expected.interface);
}
#[test]
fn g2_profile_modes_feed_back_through_the_existing_cavity_energy_owner() {
    let g=CancelGate::new_clock_free();let mut b=budget(6);b.maximum_frequency_hz=0.55;b.maximum_modes=8;
    let air=cavity(cylinder(6),b);
    let basis=air.sample(&[[0.5,0.,1.]],8,&g).unwrap();let count=basis.omegas.len();
    let mut overlap=vec![0.;count];for (c,psi) in overlap.iter_mut().zip(&basis.interface) {*c=0.1*psi[0];}
    let (system,probe)=CavityCoupling::new(&basis,1,&overlap,&vec![0.;count]).unwrap().build(
        vec![ImpactBody {potential:BodyPotential::Linear(vec![4.0]),damping_per_s:vec![0.0],
            initial:vec![ModalAcousticState {displacement_m_sqrt_kg:0.001,velocity_m_sqrt_kg_per_s:0.0}]}],
        vec![],vec![],ImpactConfig {dt_s:0.001,max_steps:100,maximum_energy_j:1.,
            energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:100.},&g).unwrap();
    let mut system=system.prepare().unwrap();let initial=system.stored_energy_j();
    let forces=vec![0.;probe.total_modes()];
    for _ in 0..100 {system.step(&forces,&g).unwrap();}
    assert!(system.state()[2..].iter().any(|v|v.abs()>1e-10));
    assert!((system.stored_energy_j()-initial).abs()<1e-10);
    let pressure=probe.pressure_at(system.state(),&air.values_at([0.5,0.,1.]).unwrap()).unwrap();
    assert!(pressure.is_finite() && pressure!=0.0);
}
