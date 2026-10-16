use super::*;
use std::f64::consts::{PI, TAU};
type Triangle = [[f64; 3]; 3];
fn box_triangles(origin: [f64; 3], size: [f64; 3]) -> Vec<Triangle> {
    let points = [[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
        [0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]]
        .map(|p| std::array::from_fn(|c| origin[c] + size[c]*p[c]));
    [[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
        [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]
        .map(|t| t.map(|i| points[i])).to_vec()
}
#[test]
fn box_spectrum_bounds_are_conservative_and_high_frequencies_keep_burton_miller() {
    let lengths = [2.,3.,4.];
    let surface = SpherePanels::from_triangles(box_triangles([0.;3], lengths)).unwrap();
    let policy = GeometryPolicy::new(&surface).unwrap();
    let exact = PI * lengths.iter().map(|l| 1./(l*l)).sum::<f64>().sqrt();
    let bound = policy.first_dirichlet_bounds()[0];
    assert!(bound <= exact && (bound-exact).abs() < 1e-13*exact);
    assert!(policy.plain_cbie_limit() < 0.8*exact);
    assert_eq!(policy.formulation(0.79*exact).unwrap(), Formulation::PlainCbie);
    assert_eq!(policy.formulation(0.81*exact).unwrap(), Formulation::BurtonMiller);
    assert_eq!(policy.formulation(exact).unwrap(), Formulation::BurtonMiller);
    for k in [0., -1., f64::NAN, f64::INFINITY] { assert!(policy.formulation(k).is_err()); }
}
#[test]
fn separation_does_not_turn_empty_space_into_a_fictitious_cavity() {
    let first = box_triangles([0.;3], [0.1,0.1,0.003]);
    let mut adjacent = first.clone(); adjacent.extend(box_triangles([0.,0.,0.1], [0.1,0.1,0.004]));
    let mut distant = first; distant.extend(box_triangles([0.,0.,10.], [0.1,0.1,0.004]));
    let a = SpherePanels::from_triangles(adjacent).unwrap(); let b = SpherePanels::from_triangles(distant).unwrap();
    let a = GeometryPolicy::new(&a).unwrap(); let b = GeometryPolicy::new(&b).unwrap();
    assert_eq!(a.first_dirichlet_bounds().len(), 2);
    assert!((a.plain_cbie_limit()-b.plain_cbie_limit()).abs() < 1e-10*a.plain_cbie_limit());
    assert_eq!(a.formulation(TAU*1000./340.).unwrap(), Formulation::PlainCbie);
    assert_eq!(b.formulation(TAU*1000./340.).unwrap(), Formulation::PlainCbie);
}
#[test]
fn axis_changes_scaling_and_malformed_solids_do_not_overstate_the_bound() {
    let raw = box_triangles([0.;3], [0.1,0.2,0.003]);
    let original = SpherePanels::from_triangles(raw.clone()).unwrap();
    let moved = SpherePanels::from_triangles(raw.iter().map(|t| t.map(|p|[3.+2.*p[1],4.+2.*p[2],5.+2.*p[0]])).collect()).unwrap();
    let a = GeometryPolicy::new(&original).unwrap(); let b = GeometryPolicy::new(&moved).unwrap();
    assert!((2.*b.plain_cbie_limit()-a.plain_cbie_limit()).abs() < 1e-11*a.plain_cbie_limit());
    let open = SpherePanels::from_triangles(raw[..raw.len()-1].to_vec()).unwrap();
    assert!(GeometryPolicy::new(&open).is_err());
    let inward = SpherePanels::from_triangles(raw.iter().map(|t|[t[0],t[2],t[1]]).collect()).unwrap();
    assert!(GeometryPolicy::new(&inward).is_err());
}
// Same dimensions and hinge pose as the previously failing native-skin/lid
// consumer. Prescribed motion is a kinematic coupon, not a piano eigenpair.
fn thin_board_and_lid() -> (SpherePanels, Vec<C64>) {
    let xy = [[0.,0.],[0.1,0.],[0.1,0.1],[0.,0.1],[0.05,0.05]];
    let top = xy.map(|p|[p[0],p[1],0.0015]); let bottom = xy.map(|p|[p[0],p[1],-0.0015]);
    let faces = [[0,1,4],[1,2,4],[2,3,4],[3,0,4]];
    let mut tris: Vec<_> = faces.iter().map(|t|t.map(|i|top[i])).collect();
    tris.extend(faces.iter().map(|t|[bottom[t[0]],bottom[t[2]],bottom[t[1]]]));
    for [a,b] in [[0,1],[1,2],[2,3],[3,0]] {
        tris.push([bottom[a],bottom[b],top[b]]); tris.push([bottom[a],top[b],top[a]]);
    }
    let lid = box_triangles([0.02,0.02,0.1], [0.05,0.05,0.004]);
    let angle = PI/6.; let (s,c) = (fs_math::det::sin(angle),fs_math::det::cos(angle));
    tris.extend(lid.iter().map(|t|t.map(|p| {
        let y=p[1]-0.02; let z=p[2]-0.1;
        [p[0],0.02+c*y-s*z,0.1+s*y+c*z]
    })));
    let mut v=vec![C64::ZERO;tris.len()]; v[..4].fill(C64::from_re(1./3.)); v[4..8].fill(C64::from_re(-1./3.));
    (SpherePanels::from_triangles(tris).unwrap(),v)
}
#[test]
fn thin_board_and_posed_lid_retain_positive_radiation_across_the_old_switch() {
    let (surface,velocity)=thin_board_and_lid(); let policy=GeometryPolicy::new(&surface).unwrap();
    let medium=Medium {density:1.2,sound_speed:340.};
    for hz in [40.,200.,280.,287.,300.] {
        let k=TAU*hz/medium.sound_speed;
        assert_eq!(policy.formulation(k).unwrap(),Formulation::PlainCbie);
        let selected=policy.solve_batch(k,medium,&[velocity.as_slice()]).unwrap().remove(0);
        let direct=helmholtz::solve_radiation(&surface,k,medium,&velocity,Formulation::PlainCbie).unwrap();
        assert_eq!(selected.pressure,direct.pressure);
        assert!(selected.radiated_power_roundoff_interval.0>0.,"at {hz} Hz: {:?}",selected.power_diagnostics);
        assert_eq!(selected.velocity,velocity);
    }
}
#[test]
fn policy_batch_keeps_phase_linearity_and_existing_refusals() {
    let (surface,velocity)=thin_board_and_lid(); let policy=GeometryPolicy::new(&surface).unwrap();
    let scale=C64::new(0.3,-0.7); let second:Vec<_>=velocity.iter().map(|v|*v*scale).collect();
    let medium=Medium {density:1.2,sound_speed:340.};let k=TAU*287./340.;
    let result=policy.solve_batch(k,medium,&[velocity.as_slice(),second.as_slice()]).unwrap();
    for (a,b) in result[0].pressure.iter().zip(&result[1].pressure) {
        assert!((*a*scale-*b).abs()<1e-11*(1.+a.abs()));
    }
    assert!(policy.solve_batch(k,medium,&[]).is_err());
    assert!(policy.solve_batch(k,medium,&[&[]]).is_err());
    assert!(policy.solve_batch(1e6,medium,&[velocity.as_slice()]).is_err());
}
