//! Source-coordinate derivatives independently checked against primal samples.
mod support;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_conduction::bc::history::march::moving::sensitivity::MovingHeatDirection;
use fs_conduction::{ConductionMesh, COMPONENT_POWER_DIMS};
use fs_conduction::fixtures::box_grid;
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cancelled_cx, with_cx};
const LENGTH: Dims = Dims([1,0,0,0,0,0]);
fn constant(value: f64, dims: Dims) -> TimeSignal { TimeSignal::Constant(QtyAny::new(value,dims)) }
fn mesh() -> ConductionMesh {
    let (complex, points) = box_grid([4,3,2],[1.0;3]);
    ConductionMesh::new(complex, points).unwrap()
}
fn source(center: [f64;3], width: [f64;3], power: f64) -> MovingHeatSource {
    MovingHeatSource::new(center.map(|v|constant(v,LENGTH)),width,constant(power,COMPONENT_POWER_DIMS),0.0,
        MovingHeatBudget{max_vertices:1000,max_elements:1000,min_active_vertices:2,power_relative_tolerance:1e-12}).unwrap()
}
fn assert_close(a: f64,b: f64,tolerance: f64) {
    assert!((a-b).abs() <= tolerance*(1.0+a.abs().max(b.abs())),"{a:e} != {b:e}");
}
#[test]
fn g3_center_width_and_power_pullbacks_match_primal_differences_on_clipped_footprints() {
    let mesh=mesh();
    with_cx(|cx| {
        for center in [[0.13,0.42,0.53],[0.55,0.68,0.30],[-0.12,0.42,0.53]] {
            let widths=[0.8,0.65,0.9];let power=3.2;
            let src=source(center,widths,power);
            let linear=src.linearize_at(cx,&mesh,0.0).unwrap();
            assert_eq!(linear.sample(),&src.at(cx,&mesh,0.0).unwrap());
            let seed:Vec<_>=mesh.positions().iter().map(|p|1.0+p[0]+2.0*p[1]*p[2]).collect();
            let gradient=linear.pullback(cx,&seed).unwrap();
            let objective=|c,w,p| {
                let sampled=source(c,w,p).at(cx,&mesh,0.0).unwrap();
                seed.iter().enumerate().map(|(i,s)|s*sampled.source.at(i)).sum::<f64>()
            };
            for control in 0..7 {
                let eps=1e-6;
                let(mut cp,mut cm,mut wp,mut wm)=(center,center,widths,widths);
                let(mut pp,mut pm)=(power,power);
                let analytic=if control<3 {cp[control]+=eps;cm[control]-=eps;gradient.center_m[control]}
                    else if control<6 {wp[control-3]+=eps;wm[control-3]-=eps;gradient.half_width_m[control-3]}
                    else {pp+=eps;pm-=eps;gradient.absorbed_power_w};
                assert_close(analytic,(objective(cp,wp,pp)-objective(cm,wm,pm))/(2.0*eps),2e-6);
            }
        }
    });
}
#[test]
fn g1_tangent_transpose_and_total_power_include_the_normalization_derivative() {
    let mesh=mesh();
    with_cx(|cx| {
        let src=source([0.13,0.42,0.53],[0.8,0.65,0.9],3.2);
        let linear=src.linearize_at(cx,&mesh,0.0).unwrap();
        let direction=MovingHeatDirection{center_m:[0.3,-0.2,0.1],half_width_m:[-0.1,0.3,0.2],absorbed_power_w:0.7};
        let applied=linear.apply(cx,direction).unwrap();
        let seed:Vec<_>=(0..mesh.vertex_count()).map(|i|(i%7)as f64-2.0).collect();
        let gradient=linear.pullback(cx,&seed).unwrap();
        let lhs:f64=seed.iter().zip(&applied).map(|(s,v)|s*v).sum();
        let rhs=gradient.absorbed_power_w*direction.absorbed_power_w
            +(0..3).map(|a|gradient.center_m[a]*direction.center_m[a]+gradient.half_width_m[a]*direction.half_width_m[a]).sum::<f64>();
        assert_close(lhs,rhs,3e-13);
        // Actual consistent P1 source-load integration, not a sum of densities.
        let mut total=0.0;let mut weights=vec![0.0;mesh.vertex_count()];
        for(e,tet)in mesh.complex().tets.iter().enumerate(){for &a in tet{for &b in tet{
            let entry=mesh.element_volume(e)/20.0*if a==b{2.0}else{1.0};
            total+=entry*applied[b as usize];weights[b as usize]+=entry;
        }}}
        assert_close(total,direction.absorbed_power_w,2e-12);
        let watts=linear.pullback(cx,&weights).unwrap();
        assert_close(watts.absorbed_power_w,1.0,2e-12);
        for axis in 0..3{assert_close(watts.center_m[axis],0.0,2e-12);assert_close(watts.half_width_m[axis],0.0,2e-12);}
    });
}
#[test]
fn g0_zero_power_derivative_requires_resolved_support_and_preserves_event_sides() {
    let mesh=mesh();
    with_cx(|cx| {
        let mut src=source([0.4,0.4,0.5],[0.8;3],0.0);
        let off=src.linearize_at(cx,&mesh,0.0).unwrap();
        let jvp=off.apply(cx,MovingHeatDirection{absorbed_power_w:1.0,..Default::default()}).unwrap();
        let unit=source([0.4,0.4,0.5],[0.8;3],1.0).at(cx,&mesh,0.0).unwrap();
        for(i,&v)in jvp.iter().enumerate(){assert_eq!(v,unit.source.at(i));}
        let unresolved=source([10.0;3],[0.8;3],0.0);
        assert!(unresolved.at(cx,&mesh,0.0).is_ok());assert!(unresolved.linearize_at(cx,&mesh,0.0).is_err());
        src=MovingHeatSource::new([constant(0.4,LENGTH),constant(0.4,LENGTH),constant(0.5,LENGTH)],[0.8;3],
            TimeSignal::Table{times:vec![0.0,1.0],values:vec![3.0,0.0],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold},
            0.0,src.budget()).unwrap();
        assert_eq!(src.linearize_at(cx,&mesh,1.0).unwrap().sample().absorbed_power_w,0.0);
        assert_eq!(src.linearize_before(cx,&mesh,1.0).unwrap().sample().absorbed_power_w,3.0);
    });
}
#[test]
fn g0_bad_seeds_directions_and_cancellation_publish_no_partial_derivative() {
    let mesh=mesh();let src=source([0.4,0.4,0.5],[0.8;3],3.0);
    let linear=with_cx(|cx|src.linearize_at(cx,&mesh,0.0).unwrap());
    with_cx(|cx|{
        assert!(linear.pullback(cx,&[]).is_err());
        assert!(linear.pullback(cx,&vec![f64::NAN;mesh.vertex_count()]).is_err());
        assert!(linear.apply(cx,MovingHeatDirection{absorbed_power_w:f64::INFINITY,..Default::default()}).is_err());
    });
    with_cancelled_cx(|cx|assert!(linear.pullback(cx,&vec![1.0;mesh.vertex_count()]).is_err()));
    with_cancelled_cx(|cx|assert!(src.linearize_at(cx,&mesh,0.0).is_err()));
    with_cx(|cx|assert!(linear.pullback(cx,&vec![1.0;mesh.vertex_count()]).is_ok()));
}
