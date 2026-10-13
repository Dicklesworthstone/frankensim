use super::*;
use crate::bridge::{Support, tests::{birth, material, temperatures, config}};

fn sweep(start: [f64;2],end: [f64;2]) -> Sweep {
    Sweep { start_m: start,end_m: end,required_clearance_m: 0.0,
        profile: Profile { flat_radius_m: 0.0002,outer_radius_m: 0.001,flank_slope: 0.5 } }
}
fn near(a: f64,b: f64) { assert!((a-b).abs()<1e-11,"{a:e} != {b:e}"); }

#[test]
fn flat_surface_clearance_margin_penetration_and_no_overlap_are_distinct() {
    let law = material(1e8,vec![]); let beam = birth(4,Support::SimplySupported,&law);
    let query = sweep([0.01,0.00051],[0.01,0.00051]);
    let clear = beam.nozzle_clearance(query,4,&mut || false).unwrap();
    near(clear.minimum.unwrap().gap_m,1e-5); assert!(!clear.interferes);
    assert!(beam.nozzle_clearance(Sweep { required_clearance_m: 2e-5,..query },4,&mut || false).unwrap().interferes);
    let inside = beam.nozzle_clearance(sweep([0.01,0.00049],[0.01,0.00049]),4,&mut || false).unwrap();
    near(inside.minimum.unwrap().gap_m,-1e-5); assert!(inside.interferes);
    let outside = beam.nozzle_clearance(sweep([-0.01,0.0],[-0.005,0.0]),4,&mut || false).unwrap();
    assert_eq!(outside.minimum,None); assert!(!outside.interferes);
}

#[test]
fn sweep_finds_interference_between_clear_motion_endpoints() {
    let law = material(1e8,vec![]); let beam = birth(4,Support::SimplySupported,&law);
    let a = [-0.01,0.00049]; let b = [0.03,0.00049];
    for p in [a,b] { assert!(!beam.nozzle_clearance(sweep(p,p),4,&mut || false).unwrap().interferes); }
    let hit = beam.nozzle_clearance(sweep(a,b),4,&mut || false).unwrap();
    assert!(hit.interferes); let w = hit.minimum.unwrap();
    assert!(w.sweep_fraction>0.0 && w.sweep_fraction<1.0); near(w.gap_m,-1e-5);
}

#[test]
fn interior_warp_and_finite_nozzle_flank_cannot_hide_behind_clear_tip() {
    let law = material(1e8,vec![]); let flat = birth(1,Support::SimplySupported,&law);
    let (beam,r) = flat.advance_controlled(1.0,&temperatures(&flat,300.0),0.0,
        &[-0.2],config(),&mut || false).unwrap();
    let tip = beam.deformed_position_m(0.004,0.0005).unwrap();
    let q = Sweep { profile: Profile { flat_radius_m: 0.0002,outer_radius_m: 0.008,flank_slope: 0.0005 },
        ..sweep([tip[0],tip[1]+2e-6],[tip[0],tip[1]+2e-6]) };
    let hit = beam.nozzle_clearance(q,1,&mut || false).unwrap();
    assert!(hit.interferes && hit.minimum.unwrap().world_x_m>tip[0]+0.001);
    let centred = Sweep { profile: Profile { flat_radius_m: 0.02,outer_radius_m: 0.02,flank_slope: 0.0 },
        ..sweep([0.01,0.0005],[0.01,0.0005]) };
    let minimum = beam.nozzle_clearance(centred,1,&mut || false).unwrap().minimum.unwrap();
    near(minimum.reference_x_m,0.01);
    near(minimum.gap_m,-r.maximum_absolute_deflection_m);
}

#[test]
fn moving_vertical_clearance_and_thermal_axial_motion_use_actual_shape() {
    let law = material(1e8,vec![]); let flat = birth(4,Support::SimplySupported,&law);
    let (beam,_) = flat.advance_controlled(1.0,&temperatures(&flat,500.0),0.0,
        &[0.0;4],config(),&mut || false).unwrap();
    let q = sweep([0.0203,0.0008],[0.0203,0.0001]);
    let report = beam.nozzle_clearance(q,4,&mut || false).unwrap();
    assert!(report.interferes);
    near(report.minimum.unwrap().sweep_fraction,1.0);
    let reversed = beam.nozzle_clearance(Sweep { start_m:q.end_m,end_m:q.start_m,..q },4,&mut || false).unwrap();
    near(report.minimum.unwrap().gap_m,reversed.minimum.unwrap().gap_m);
    let end = beam.deformed_position_m(beam.length_m(),0.0005).unwrap()[0];
    near(report.minimum.unwrap().world_x_m,end);
    // The expanded surface reaches farther into the nozzle's rising flank.
    let prior = flat.nozzle_clearance(q,4,&mut || false).unwrap();
    assert!(report.minimum.unwrap().gap_m < prior.minimum.unwrap().gap_m);
}

#[test]
fn cancelled_incomplete_or_invalid_queries_return_no_clearance_claim() {
    let law = material(1e8,vec![]); let beam = birth(4,Support::SimplySupported,&law);
    let q = sweep([0.0,0.001],[0.02,0.001]); let before = beam.clone();
    assert!(beam.nozzle_clearance(q,3,&mut || false).is_err());
    let mut n=0;
    assert!(beam.nozzle_clearance(q,4,&mut || { n+=1; n>3 }).is_err());
    for p in [Profile { flat_radius_m:0.002,..q.profile },Profile { flank_slope:-1.0,..q.profile },
        Profile { outer_radius_m:f64::NAN,..q.profile }] {
        assert!(beam.nozzle_clearance(Sweep { profile:p,..q },4,&mut || false).is_err());
    }
    assert_eq!(beam,before);
}
