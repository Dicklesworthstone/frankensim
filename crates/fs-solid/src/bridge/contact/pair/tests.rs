use super::*;
use crate::bridge::{Support, tests::{birth, config, material, temperatures}};

fn query() -> Query {
    Query { lower_origin_m: [0.0,0.0], upper_origin_m: [0.0,0.001],
        width_m: 0.001, capture_gap_m: 1e-10, penetration_tolerance_m: 1e-11 }
}
fn budget() -> Budget { Budget { max_cells: 256, max_rows: 256 } }
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() < tolerance, "{a:e} != {b:e}");
}
fn totals(f: &Footprint) -> [f64; 3] {
    f.rows.iter().fold([0.0;3], |mut sum,r| {
        sum[0] += r.lower_reference_area_m2;
        sum[1] += r.upper_reference_area_m2;
        sum[2] += r.projected_area_m2;
        sum
    })
}

// G0/G3: different reference partitions meet in a shared world partition.
#[test]
fn translated_nonmatching_meshes_preserve_both_material_areas() {
    let law = material(1e8,vec![]);
    let lower = birth(3,Support::SimplySupported,&law);
    let upper = birth(5,Support::SimplySupported,&law);
    let q = Query { upper_origin_m: [0.001,0.001], ..query() };
    let result = lower.neighbour_contact(&upper,q,budget(),&mut || false).unwrap();
    close(result.overlap_projected_length_m,0.019,1e-14);
    assert_eq!(result.cells_visited,8);
    assert!(result.rows.len() <= 7);
    for value in totals(&result) { close(value,0.019*0.001,1e-16); }
    for row in &result.rows {
        assert!(row.world_interval_m[1] > row.world_interval_m[0]);
        close(row.lower_reference_interval_m[0],row.world_interval_m[0],1e-14);
        close(row.upper_reference_interval_m[0]+0.001,row.world_interval_m[0],1e-14);
    }
    for pair in result.rows.windows(2) {
        close(pair[0].world_interval_m[1],pair[1].world_interval_m[0],1e-14);
    }
    let moved = lower.neighbour_contact(&upper,Query { lower_origin_m: [0.1,0.1],
        upper_origin_m: [0.101,0.101], ..q },budget(),&mut || false).unwrap();
    for (a,b) in totals(&result).into_iter().zip(totals(&moved)) { close(a,b,1e-15); }
}

// G1: two affine surface stretches have distinct reference-area Jacobians.
#[test]
fn differential_thermal_expansion_never_becomes_one_shared_reference_area() {
    let law = material(1e8,vec![]);
    let lower = birth(7,Support::SimplySupported,&law);
    let upper = birth(11,Support::SimplySupported,&law);
    let lower = lower.advance_controlled(1.0,&temperatures(&lower,500.0),0.0,
        &[0.0;7],config(),&mut || false).unwrap().0;
    let upper = upper.advance_controlled(1.0,&temperatures(&upper,300.0),0.0,
        &[0.0;11],config(),&mut || false).unwrap().0;
    let result = lower.neighbour_contact(&upper,query(),budget(),&mut || false).unwrap();
    let areas = totals(&result);
    close(areas[0],0.001*0.02/1.002,1e-14);
    close(areas[1],0.001*0.02,1e-14);
    close(areas[2],areas[1],1e-15);
    assert!(areas[0] < areas[1]);
    assert_eq!(lower.time_s(),1.0);
    assert_eq!(upper.time_s(),1.0);
}

// G1: the same quadratic capture limit as a flat later pass, but here the upper
// surface is an actual separate beam with a different discretization.
#[test]
fn relative_curvature_resolves_internal_gap_and_capture_roots() {
    let law = material(1e8,vec![]);
    let lower = birth(1,Support::SimplySupported,&law);
    let upper = birth(5,Support::SimplySupported,&law);
    let (lower,r) = lower.advance_controlled(1.0,&temperatures(&lower,300.0),9.81,
        &[0.0],config(),&mut || false).unwrap();
    let upper = upper.advance_controlled(1.0,&temperatures(&upper,300.0),0.0,
        &[0.0;5],config(),&mut || false).unwrap().0;
    let q = Query { capture_gap_m: r.maximum_sag_m/4.0,..query() };
    let result = lower.neighbour_contact(&upper,q,budget(),&mut || false).unwrap();
    let fraction = 1.0-(0.75_f64).sqrt();
    close(totals(&result)[0],0.001*0.02*fraction,1e-14);
    assert!(totals(&result)[1] < totals(&result)[0]);
    close(result.rows.iter().map(|r| r.maximum_gap_m).fold(0.0,f64::max),r.maximum_sag_m,1e-13);
    assert!(result.rows.iter().any(|r| r.projected_area_m2 == 0.0));
}

#[test]
fn interior_penetration_is_not_hidden_by_clear_boundary_points() {
    let law = material(1e8,vec![]);
    let lower = birth(3,Support::SimplySupported,&law);
    let upper = birth(1,Support::SimplySupported,&law);
    let lower = lower.advance_controlled(1.0,&temperatures(&lower,300.0),0.0,
        &[0.0;3],config(),&mut || false).unwrap().0;
    let (upper,r) = upper.advance_controlled(1.0,&temperatures(&upper,300.0),9.81,
        &[0.0],config(),&mut || false).unwrap();
    let before = (lower.clone(),upper.clone());
    assert!(matches!(lower.neighbour_contact(&upper,query(),budget(),&mut || false),
        Err(Error::Invalid("deformed neighbour faces penetrate between stations"))));
    let clear = lower.neighbour_contact(&upper,Query {
        upper_origin_m: [0.0,0.001+2.0*r.maximum_sag_m],..query()
    },budget(),&mut || false).unwrap();
    assert!(clear.minimum_gap_m.unwrap() > 0.0);
    assert_eq!(totals(&clear),[0.0;3]);
    assert_eq!((lower,upper),before);
}

#[test]
fn no_overlap_is_distinct_from_zero_captured_area_and_still_budgeted() {
    let law = material(1e8,vec![]);
    let lower = birth(3,Support::SimplySupported,&law);
    let upper = birth(5,Support::SimplySupported,&law);
    let q = Query { upper_origin_m: [0.1,0.001],..query() };
    let result = lower.neighbour_contact(&upper,q,budget(),&mut || false).unwrap();
    assert!(result.rows.is_empty());
    assert_eq!(result.minimum_gap_m,None);
    assert_eq!(result.overlap_projected_length_m,0.0);
    assert_eq!(result.cells_visited,8);
    assert!(lower.neighbour_contact(&upper,q,Budget { max_cells:7,..budget() },&mut || false).is_err());
    assert!(lower.neighbour_contact(&upper,query(),Budget { max_rows:0,..budget() },&mut || false).is_err());
}

#[test]
fn time_width_nonfinite_and_late_cancellation_refuse_without_partial_rows() {
    let law = material(1e8,vec![]);
    let lower = birth(3,Support::SimplySupported,&law);
    let upper = birth(5,Support::SimplySupported,&law);
    for q in [Query { width_m:0.002,..query() }, Query { capture_gap_m:f64::NAN,..query() },
        Query { upper_origin_m:[f64::INFINITY,0.001],..query() }] {
        assert!(lower.neighbour_contact(&upper,q,budget(),&mut || false).is_err());
    }
    let stale = upper.advance_controlled(1.0,&temperatures(&upper,300.0),0.0,
        &[0.0;5],config(),&mut || false).unwrap().0;
    assert!(lower.neighbour_contact(&stale,query(),budget(),&mut || false).is_err());
    let mut polls = 0;
    lower.neighbour_contact(&upper,query(),budget(),&mut || { polls+=1; false }).unwrap();
    let mut calls = 0;
    assert_eq!(lower.neighbour_contact(&upper,query(),budget(),&mut || {
        calls+=1; calls>=polls
    }),Err(Error::Cancelled));
    assert_eq!(lower.time_s(),0.0);
    assert_eq!(upper.time_s(),0.0);
}
