use super::*;
use crate::helmholtz::{Formulation, solve_radiation};

fn box_triangles(origin:Point,size:Point) -> Vec<Triangle> {
    let points=[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
        [0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]]
        .map(|p|std::array::from_fn(|i|origin[i]+size[i]*p[i]));
    [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
        [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]
        .map(|t|t.map(|i|points[i])).to_vec()
}
fn slab() -> SpherePanels {
    SpherePanels::from_triangles(box_triangles([0.,0.,-0.002],[1.,1.,0.004])).unwrap()
}
#[test]
fn close_exterior_points_inside_the_bounding_sphere_are_admitted_by_actual_surface_distance() {
    let body=slab();let points=[[0.3,0.4,0.022],[0.3,0.4,-0.022]];
    let geometry=Geometry::new(&body,&points,0.01).unwrap();
    for &d in geometry.clearances_m() {assert!((d-0.02).abs()<1e-14);}
    assert_eq!(geometry.points_m(),&points);
    for point in [[0.5,0.5,0.],[0.5,0.5,0.002],[0.,0.,0.002],
        [0.,0.5,0.002],[0.5,0.5,0.002001],[f64::NAN,0.,0.]] {
        assert!(Geometry::new(&body,&[point],0.01).is_err(),"accepted {point:?}");
    }
    let mut pair=box_triangles([0.,0.,0.],[0.1;3]);pair.extend(box_triangles([0.3,0.,0.],[0.1;3]));
    let pair=SpherePanels::from_triangles(pair).unwrap();
    assert!(Geometry::new(&pair,&[[0.2,0.05,0.05]],0.01).is_ok());
    assert!(Geometry::new(&pair,&[[0.35,0.05,0.05]],0.01).is_err());
}
#[test]
fn closed_triangle_admission_and_clearance_are_rigid_frame_covariant() {
    let original=box_triangles([0.,0.,0.],[0.1,0.2,0.003]);
    let p=[0.03,0.04,0.023];
    let angle=0.713_f64;let shift=[12.,-8.,3.];
    let pose=|p:Point|[shift[0]+angle.cos()*p[0]-angle.sin()*p[2],shift[1]+p[1],
        shift[2]+angle.sin()*p[0]+angle.cos()*p[2]];
    let body=SpherePanels::from_triangles(original.iter().map(|t|t.map(pose)).collect()).unwrap();
    assert!((Geometry::new(&body,&[pose(p)],0.01).unwrap().clearances_m()[0]-0.02).abs()<1e-12);
    let mut open=original.clone();open.pop();
    assert!(Geometry::new(&SpherePanels::from_triangles(open).unwrap(),&[p],0.01).is_err());
    let reversed=original.into_iter().map(|[a,b,c]|[a,c,b]).collect();
    assert!(Geometry::new(&SpherePanels::from_triangles(reversed).unwrap(),&[p],0.01).is_err());
    let legacy=SpherePanels::new(body.centroids().to_vec(),body.normals().to_vec(),body.areas().to_vec()).unwrap();
    assert!(Geometry::new(&legacy,&[pose(p)],0.01).is_err());
}
#[test]
fn near_triangle_static_double_layer_matches_the_exact_solid_angle_not_a_centroid_proxy() {
    let t=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]];
    for x in [[0.3,0.3,0.01],[0.5,0.,0.02],[-0.01,-0.01,0.03],[0.2,0.3,-0.015]] {
        let mut work=Work {options:Options::default(),evaluations:0};
        let integral=integrate(0.,x,t,[0.,0.,1.],0,&mut work).unwrap();
        let exact=-angle(x,t)/FOUR_PI;
        assert!((integral.d.re-exact).abs()<2e-8,"{} vs {exact}",integral.d.re);
        assert_eq!(integral.d.im,0.);assert!(integral.s.re>0.);
        assert!(work.evaluations>80);
    }
    // Uniform double layer has exactly zero exterior potential on a closed
    // body, including close points that fooled centroid winding evaluation.
    let body=slab();let x=[0.3,0.4,0.022];let mut total=C64::ZERO;
    let mut work=Work {options:Options::default(),evaluations:0};
    for (&t,&n) in body.triangles().unwrap().iter().zip(body.normals()) {
        total=total+integrate(0.,x,t,n,0,&mut work).unwrap().d;
    }
    assert!(total.abs()<5e-8,"{total:?}");
}
#[test]
fn prepared_rows_observe_the_same_bem_fields_with_linearity_far_limit_and_strict_cross_wiring() {
    let surface=SpherePanels::from_triangles(box_triangles([0.,0.,-0.002],[0.1,0.1,0.004])).unwrap();
    let medium=Medium {density:1.2,sound_speed:340.};let k=std::f64::consts::TAU*100./340.;
    let velocity:Vec<_>=surface.normals().iter().map(|n|C64::from_re(n[2])).collect();
    let solution=solve_radiation(&surface,k,medium,&velocity,Formulation::PlainCbie).unwrap();
    let points=[[0.04,0.04,0.022],[0.04,0.04,2.]];
    let geometry=Geometry::new(&surface,&points,0.01).unwrap();
    let prepared=geometry.prepare(k,medium,Options::default()).unwrap();
    let observation=prepared.evaluate(&solution).unwrap();
    assert!(observation.pressure.iter().all(|p|p.abs()>1e-12));
    assert!(observation.quadrature_error_estimate_pa.iter().all(|e|e.is_finite()&&*e>=0.));
    let far=helmholtz::exterior_pressure_at_points(&surface,&solution,medium,&points[1..]).unwrap();
    assert!((observation.pressure[1]-far[0]).abs()<0.01*far[0].abs());
    let gain=C64::new(0.3,-0.4);let mut scaled=solution.clone();
    for value in scaled.pressure.iter_mut().chain(&mut scaled.velocity) {*value=*value*gain;}
    let changed=prepared.evaluate(&scaled).unwrap();
    for (&a,&b) in observation.pressure.iter().zip(&changed.pressure) {assert!((b-a*gain).abs()<1e-12*(1.+a.abs()));}
    assert_eq!(prepared.evaluate(&solution).unwrap().pressure,observation.pressure);
    scaled=solution.clone();scaled.k*=2.;assert!(prepared.evaluate(&scaled).is_err());
    scaled=solution.clone();scaled.medium.density*=2.;assert!(prepared.evaluate(&scaled).is_err());
    scaled=solution.clone();scaled.surface_fingerprint^=1;assert!(prepared.evaluate(&scaled).is_err());
    scaled=solution;scaled.pressure[0]=C64::new(f64::NAN,0.);assert!(prepared.evaluate(&scaled).is_err());
}
#[test]
fn under_resolved_work_refuses_instead_of_publishing_centroid_or_partial_receiver_results() {
    let body=slab();let geometry=Geometry::new(&body,&[[0.3,0.4,0.022],[2.,2.,2.]],0.01).unwrap();
    let medium=Medium {density:1.2,sound_speed:340.};
    assert!(geometry.prepare(2.,medium,Options {maximum_kernel_evaluations:80,..Options::default()}).is_err());
    assert!(geometry.prepare(2.,medium,Options {maximum_depth:0,..Options::default()}).is_err());
    assert!(geometry.prepare(2.,medium,Options {relative_tolerance:f64::NAN,..Options::default()}).is_err());
    assert!(geometry.prepare(f64::NAN,medium,Options::default()).is_err());
    assert_eq!(geometry.points_m(),&[[0.3,0.4,0.022],[2.,2.,2.]]);
}
