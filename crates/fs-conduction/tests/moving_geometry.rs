//! G0/G3: authored geometry edits preserve the original source physics.
mod support;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::moving::sensitivity::geometry::MovingHeatGeometryControl as G;
use fs_conduction::{ConductionMesh, COMPONENT_POWER_DIMS};
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::{TimeSignal, Interp};
use support::{with_cx, with_cancelled_cx};

fn source() -> MovingHeatSource {
    let length = Dims([1,0,0,0,0,0]);
    MovingHeatSource::new([
        TimeSignal::Ramp { t_start:0., t_end:2., from:QtyAny::new(0.1,length), to:QtyAny::new(0.7,length) },
        TimeSignal::Table { times:vec![0.,1.,2.], values:vec![0.2,0.5,0.8], dims:length, interp:Interp::Hold },
        TimeSignal::Constant(QtyAny::new(0.3,length)),
    ], [2.;3], TimeSignal::Table { times:vec![0.,1.,2.], values:vec![2.,3.,0.],
        dims:COMPONENT_POWER_DIMS, interp:Interp::Hold }, 0.,
        MovingHeatBudget { max_vertices:27,max_elements:48,min_active_vertices:1,power_relative_tolerance:1e-11 }).unwrap()
}
fn mesh() -> ConductionMesh {
    let (c,p)=fs_conduction::fixtures::box_grid([2;3],[1.;3]);
    ConductionMesh::new(c,p).unwrap()
}
#[test]
fn g3_geometry_edits_keep_knots_event_sides_power_and_the_original_source() {
    let original=source(); let saved=original.clone(); let m=mesh();
    with_cx(|cx| {
        let changed=original.with_geometry_values(cx,&[
            (G::CenterValue {axis:0,value:1},-0.1),
            (G::CenterValue {axis:1,value:1},0.4),
            (G::CenterValue {axis:2,value:0},0.),
            (G::HalfWidth {axis:0},1.5),
        ],0.,32).unwrap();
        assert_eq!(original,saved); assert_eq!(changed.budget(),original.budget());
        assert_eq!(changed.absorbed_power_values(3).unwrap(),vec![2.,3.,0.]);
        let right=changed.at(cx,&m,1.).unwrap(); let left=changed.before(cx,&m,1.).unwrap();
        assert_eq!(right.center_m,[0.,0.4,0.]); assert_eq!(left.center_m,[0.,0.2,0.]);
        assert_eq!(right.absorbed_power_w,3.); assert_eq!(left.absorbed_power_w,2.);
        for sample in [right,left] {
            assert!((sample.delivered_power_w-sample.absorbed_power_w).abs()<1e-11);
        }
        assert_eq!(changed.geometry_value(G::HalfWidth {axis:0}).unwrap(),1.5);
        assert_eq!(changed.geometry_value(G::HalfWidth {axis:1}).unwrap(),2.);
        assert_eq!(changed.at(cx,&m,2.).unwrap().absorbed_power_w,0.);
    });
}
#[test]
fn g0_malformed_batches_do_not_publish_partial_source_edits() {
    let source=source(); let original=source.clone();
    with_cx(|cx| {
        for edits in [vec![(G::CenterValue {axis:3,value:0},1.)],
            vec![(G::CenterValue {axis:0,value:2},1.)],
            vec![(G::HalfWidth {axis:1},0.)],vec![(G::HalfWidth {axis:0},f64::INFINITY)],
            vec![(G::CenterValue {axis:2,value:0},f64::NAN)],
            vec![(G::HalfWidth {axis:0},1.),(G::HalfWidth {axis:0},2.)]] {
            assert!(source.with_geometry_values(cx,&edits,0.,32).is_err());
            assert_eq!(source,original);
        }
        assert!(source.with_geometry_values(cx,&[],0.,1).is_err());
        assert!(source.with_geometry_values(cx,&[],f64::NAN,32).is_err());
        assert_eq!(source.with_geometry_values(cx,&[],0.,32).unwrap(),source);
    });
    with_cancelled_cx(|cx| assert!(source.with_geometry_values(cx,&[(G::HalfWidth {axis:0},1.)],0.,32).is_err()));
    assert_eq!(source,original);
}
