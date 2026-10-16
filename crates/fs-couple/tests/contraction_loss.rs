use fs_couple::bernoulli_aperture::contraction::{ContractionLoss, ContractionSpec};
use fs_material::gas::{GasSpec, GasState};
fn gas(pressure:f64)->GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,pressure).unwrap()
}
fn spec(ratio:f64)->ContractionSpec {
    ContractionSpec { area_ratio:ratio, maximum_contracted_mach:0.2,
        provenance:"authored symmetric separated-window comparison".into() }
}
fn loss(ratio:f64)->ContractionLoss {ContractionLoss::new(0.01,0.004,gas(101325.0),spec(ratio)).unwrap()}

#[test]
fn bernoulli_inverse_is_odd_and_has_no_artificial_pressure_dead_zone() {
    let l=loss(0.6);
    for dp in [1e-40_f64,1e-20,1e-12,1.0,100.0] {
        for sign in [-1.0,1.0] {
            let q=sign*0.6*0.01*0.004*(2.0*dp/l.gas().density).sqrt();
            let f=l.at_flow(q).unwrap();
            assert!((f.pressure_drop_pa-sign*dp).abs()<2e-12*dp);
            assert!(f.dissipated_power_w>0.0);
            let solved=l.solve(0.0,sign*dp).unwrap();
            assert!((solved.volume_flow_m3_s-q).abs()<2e-12*q.abs());
        }
    }
    assert_eq!(l.at_flow(0.0).unwrap().dissipated_power_w,0.0);
}

#[test]
fn geometry_density_and_contraction_change_only_the_original_bernoulli_scaling() {
    let a=loss(0.6); let q=1e-5;
    let base=a.at_flow(q).unwrap();
    let doubled_flow=a.at_flow(2.0*q).unwrap();
    assert!((doubled_flow.pressure_drop_pa/base.pressure_drop_pa-4.0).abs()<1e-12);
    assert!((doubled_flow.dissipated_power_w/base.dissipated_power_w-8.0).abs()<1e-12);
    let wider=ContractionLoss::new(0.02,0.004,*a.gas(),spec(0.6)).unwrap();
    let denser=ContractionLoss::new(0.01,0.004,gas(202650.0),spec(0.6)).unwrap();
    assert!((wider.at_flow(q).unwrap().pressure_drop_pa/base.pressure_drop_pa-0.25).abs()<1e-12);
    assert!((denser.at_flow(q).unwrap().pressure_drop_pa/base.pressure_drop_pa-2.0).abs()<1e-12);
    assert!((loss(0.3).at_flow(q).unwrap().pressure_drop_pa/base.pressure_drop_pa-4.0).abs()<1e-12);
    assert_eq!(a.gas(),&gas(101325.0));
    assert_eq!(a.spec(),&spec(0.6));
}

#[test]
fn simultaneous_linear_and_nonlinear_load_matches_manufactured_signed_flows() {
    let l=loss(0.6);
    for q in [-1e-4_f64,-1e-8,0.0,1e-8,1e-4] {
        for r in [0.0,1e-4,1e4,1e6,1e10] {
            let drop=0.5*l.gas().density*q*q.abs()/(0.6*0.01*0.004_f64).powi(2);
            let f=l.solve(r,r*q+drop).unwrap();
            assert!((f.volume_flow_m3_s-q).abs()<2e-12*q.abs().max(1e-30));
            assert!(f.dissipated_power_w>=0.0);
        }
    }
}

#[test]
fn bad_material_geometry_or_solved_mach_refuses_without_clipping() {
    let l=loss(0.6);
    for bad in [0.0,-1.0,1.1,f64::NAN] {
        assert!(ContractionLoss::new(0.01,0.004,gas(101325.0),spec(bad)).is_err());
    }
    let mut s=spec(0.6);s.provenance.clear();
    assert!(ContractionLoss::new(0.01,0.004,gas(101325.0),s).is_err());
    let mut wrong=gas(101325.0);wrong.density*=2.0;
    assert!(ContractionLoss::new(0.01,0.004,wrong,spec(0.6)).is_err());
    assert!(ContractionLoss::new(0.0,0.004,gas(101325.0),spec(0.6)).is_err());
    let limit=l.effective_area_m2()*l.gas().sound_speed*0.2;
    assert!(l.at_flow(limit*1.01).is_err());
    assert!(l.solve(0.0,1e6).is_err());
    assert!(l.solve(-1.0,0.0).is_err());
    assert!(l.solve(1.0,f64::INFINITY).is_err());
    assert!(l.at_flow(f64::NAN).is_err());
}
