//! G0/G1/G4/G5: oriented implicit interfaces and unchanged bulk rules.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::quad3::{cut_cell_rules3, cut_cell_rules3_with_surface, QuadratureControl3, QuadratureError3, QuadratureOptions3, QuadratureWork3};
use fs_cutfem::quad3::surface::{surface_cell_rules3, SurfaceOptions3, SurfaceRules3};
use fs_ivl::Interval;
struct Plane { n: [f64; 3], d: f64 }
fn axis(a: HeightAxis) -> usize { match a { HeightAxis::X => 0, HeightAxis::Y => 1, HeightAxis::Z => 2 } }
impl CutSdf3 for Plane {
    fn value(&self, p: [f64; 3]) -> f64 { self.n.iter().zip(p).map(|(n,x)| n*x).sum::<f64>() - self.d }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let mut v = Interval::new(-self.d, -self.d);
        for i in 0..3 { v = v + Interval::new(self.n[i], self.n[i])*Interval::new(lo[i],hi[i]); }
        v
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval { let n=self.n[axis(a)]; Interval::new(n,n) }
}
struct Sphere;
impl CutSdf3 for Sphere {
    fn value(&self,p:[f64;3])->f64 { p.iter().map(|x|(x-0.5)*(x-0.5)).sum::<f64>()-0.35*0.35 }
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let r=Interval::new(0.35,0.35);let mut value=-(r*r);
        for i in 0..3 { let d=Interval::new(lo[i],hi[i])-Interval::new(0.5,0.5); value=value+d*d; } value
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval {
        let i=axis(a); Interval::new(2.0,2.0)*(Interval::new(lo[i],hi[i])-Interval::new(0.5,0.5))
    }
}
fn unit()->HexCell { HexCell::try_new([0.0;3],[1.0;3]).unwrap() }
fn rules(sdf:&dyn CutSdf3,cell:HexCell,depth:u32)->SurfaceRules3 {
    let mut poll=|_|ControlFlow::Continue(());
    let mut c=QuadratureControl3::new(QuadratureOptions3{depth,..Default::default()},&mut poll).unwrap();
    surface_cell_rules3(sdf,cell,Default::default(),&mut c).unwrap()
}
#[test]
fn g0_planar_area_normal_moments_and_positive_weights_in_every_direction() {
    for i in 0..3 { for sign in [-1.0,1.0] {
        let mut n=[0.0;3];n[i]=sign;let p=Plane{n,d:sign*0.37};let rule=rules(&p,unit(),2);
        assert!((rule.area()-1.0).abs()<1e-12);
        for node in rule.points() {
            assert!(node.weight>0.0);assert!((node.position[i]-0.37).abs()<1e-9);assert_eq!(node.normal,n);
        }
        for a in 0..3 {let moment:f64=rule.points().iter().map(|p|p.weight*p.position[a]).sum();
            assert!((moment-if a==i{0.37}else{0.5}).abs()<1e-9);
        }
    } }
}
#[test]
fn g0_coincident_subdivision_and_background_faces_have_one_material_side_owner() {
    let left=HexCell::try_new([0.0;3],[0.5,1.0,1.0]).unwrap();
    let right=HexCell::try_new([0.5,0.0,0.0],[1.0;3]).unwrap();
    for sign in [-1.0,1.0] {for depth in 0..=3 {
        let plane=Plane{n:[sign,0.0,0.0],d:sign*0.5};
        let a=rules(&plane,left,depth);let b=rules(&plane,right,depth);
        assert!((a.area()+b.area()-1.0).abs()<1e-12);
        assert!((rules(&plane,unit(),depth).area()-1.0).abs()<1e-12);
        if sign>0.0 {assert!(b.points().is_empty());}else{assert!(a.points().is_empty());}
    }}
    assert!((rules(&Plane{n:[1.0,0.0,0.0],d:1.0},unit(),2).area()-1.0).abs()<1e-12);
    assert!(rules(&Plane{n:[-1.0,0.0,0.0],d:-1.0},unit(),2).points().is_empty());
}
#[test]
fn g1_tilted_graph_uses_surface_measure_not_projected_base_area() {
    let p=Plane{n:[-0.2,-0.1,1.0],d:0.2};let r=rules(&p,unit(),0);let exact=1.05_f64.sqrt();
    assert!((r.area()-exact).abs()<1e-12);
    for node in r.points(){for i in 0..3{assert!((node.normal[i]-p.n[i]/exact).abs()<1e-12);}}
    let moment:f64=r.points().iter().map(|p|p.weight*p.position[2]).sum();
    assert!((moment-0.35*exact).abs()<1e-9);
}
#[test]
fn g1_sphere_has_oriented_area_and_balanced_constant_pressure_resultant() {
    let r=rules(&Sphere,unit(),4);let exact=4.0*std::f64::consts::PI*0.35*0.35;
    assert!((r.area()-exact).abs()/exact<0.02,"area {} vs {exact}",r.area());
    let mut resultant=[0.0;3];let mut divergence=0.0;
    for p in r.points(){
        assert!(Sphere.value(p.position).abs()<1e-9);
        for i in 0..3{resultant[i]+=p.weight*p.normal[i];divergence+=p.weight*(p.position[i]-0.5)*p.normal[i]/3.0;}
    }
    assert!(resultant.iter().all(|v|v.abs()<1e-8));
    assert!((divergence-r.area()*0.35/3.0).abs()<1e-9);
}
#[test]
fn g3_opt_in_surface_retention_does_not_change_bulk_numbers() {
    let plane=Plane{n:[1.0,0.0,0.0],d:0.37};
    let mut cp=|_|ControlFlow::Continue(());let mut c=QuadratureControl3::new(Default::default(),&mut cp).unwrap();
    let old=cut_cell_rules3(&plane,unit(),&mut c).unwrap();let old_work=c.work();
    let mut cp=|_|ControlFlow::Continue(());let mut c=QuadratureControl3::new(Default::default(),&mut cp).unwrap();
    let new=cut_cell_rules3_with_surface(&plane,unit(),Default::default(),&mut c).unwrap();
    assert_eq!(old.bulk(),new.bulk());assert_eq!(old.volume_bounds(),new.volume_bounds());
    assert!(old.surface().is_none());assert!((new.surface().unwrap().area()-1.0).abs()<1e-12);
    assert!(c.work().points>old_work.points);assert!(c.work().field_evaluations>old_work.field_evaluations);
}
#[test]
fn g4_surface_budgets_and_cancellation_never_publish_a_partial_rule() {
    let plane=Plane{n:[1.0,0.0,0.0],d:0.37};
    let mut cp=|_|ControlFlow::Continue(());let mut c=QuadratureControl3::new(QuadratureOptions3{depth:0,max_points:9,..Default::default()},&mut cp).unwrap();
    surface_cell_rules3(&plane,unit(),Default::default(),&mut c).unwrap();
    assert!(matches!(surface_cell_rules3(&plane,unit(),Default::default(),&mut c),Err(QuadratureError3::PointBudget)));assert_eq!(c.work().points,9);
    let mut cp=|w:QuadratureWork3|if w.points>0{ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=QuadratureControl3::new(QuadratureOptions3{depth:0,..Default::default()},&mut cp).unwrap();
    assert!(matches!(surface_cell_rules3(&plane,unit(),Default::default(),&mut c),Err(QuadratureError3::Cancelled)));assert_eq!(c.work().points,1);
    let mut cp=|_|ControlFlow::Continue(());let mut c=QuadratureControl3::new(QuadratureOptions3{depth:0,..Default::default()},&mut cp).unwrap();
    assert!(matches!(surface_cell_rules3(&Sphere,unit(),Default::default(),&mut c),Err(QuadratureError3::UnresolvedCell(_))));
}
#[test]
fn g0_inconsistent_point_samples_and_uncertain_normals_are_refused() {
    struct Bad { broad:bool }
    impl CutSdf3 for Bad {
        fn value(&self,p:[f64;3])->f64{p[0]-0.37+if self.broad{0.0}else{1.0}}
        fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{Interval::new(lo[0],hi[0])-Interval::new(0.37,0.37)}
        fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval{
            if a==HeightAxis::X {if self.broad {Interval::new(0.9,1.1)}else{Interval::new(1.0,1.0)}}else{Interval::new(0.0,0.0)}
        }
    }
    for broad in [false,true]{let mut cp=|_|ControlFlow::Continue(());let mut c=QuadratureControl3::new(Default::default(),&mut cp).unwrap();
        assert!(matches!(surface_cell_rules3(&Bad{broad},unit(),Default::default(),&mut c),Err(QuadratureError3::Invalid(_))));}
}
#[test]
fn g5_replay_empty_domains_and_no_implicit_box_surface() {
    let p=Plane{n:[1.0,0.0,0.0],d:0.37};assert_eq!(rules(&p,unit(),2).points(),rules(&p,unit(),2).points());
    for d in [-1.0,2.0]{assert!(rules(&Plane{n:[1.0,0.0,0.0],d},unit(),2).points().is_empty());}
    let mut cp=|_|ControlFlow::Continue(());let mut c=QuadratureControl3::new(Default::default(),&mut cp).unwrap();
    assert!(surface_cell_rules3(&p,unit(),SurfaceOptions3{relative_normal_tolerance:f64::NAN},&mut c).is_err());assert_eq!(c.work().boxes,0);
}
