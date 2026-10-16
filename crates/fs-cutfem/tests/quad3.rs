use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::quad3::{cut_cell_rules3, CutRules3, QuadratureControl3, QuadratureError3, QuadratureOptions3};
use fs_ivl::Interval;

struct Plane { normal: [f64; 3], offset: f64 }
impl CutSdf3 for Plane {
    fn value(&self, p: [f64; 3]) -> f64 {
        p.iter().zip(self.normal).map(|(x, n)| x*n).sum::<f64>() - self.offset
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let mut v = Interval::new(-self.offset, -self.offset);
        for a in 0..3 { v = v + Interval::new(lo[a], hi[a]) * Interval::new(self.normal[a], self.normal[a]); }
        v
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let a = match a { HeightAxis::X => 0, HeightAxis::Y => 1, HeightAxis::Z => 2 };
        Interval::new(self.normal[a], self.normal[a])
    }
}
struct Sphere;
impl CutSdf3 for Sphere {
    fn value(&self, p: [f64; 3]) -> f64 { p.iter().map(|x| (x-0.5)*(x-0.5)).sum::<f64>() - 0.35*0.35 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let r = Interval::new(0.35, 0.35);
        let mut v = -(r*r);
        for a in 0..3 {
            let d = Interval::new(lo[a], hi[a]) - Interval::new(0.5, 0.5);
            v = v + d*d;
        }
        v
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        let a = match a { HeightAxis::X => 0, HeightAxis::Y => 1, HeightAxis::Z => 2 };
        Interval::new(2.0, 2.0) * (Interval::new(lo[a], hi[a]) - Interval::new(0.5, 0.5))
    }
}
fn unit() -> HexCell { HexCell::try_new([0.0; 3], [1.0; 3]).unwrap() }
fn rule(sdf: &dyn CutSdf3, depth: u32) -> CutRules3 {
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { depth, ..Default::default() }, &mut callback).unwrap();
    cut_cell_rules3(sdf, unit(), &mut control).unwrap()
}
#[test]
fn slab_volume_moments_and_normal_reversal_use_positive_weights() {
    for axis in 0..3 {
        let mut n = [0.0; 3]; n[axis] = 1.0;
        let p = Plane { normal: n, offset: 0.37 };
        let q = rule(&p, 1);
        assert!((q.volume() - 0.37).abs() < 1e-9);
        let moment: f64 = q.bulk().iter().map(|(x,w)| x[axis]*w).sum();
        assert!((moment - 0.37*0.37/2.0).abs() < 1e-9);
        assert!(q.bulk().iter().all(|(_,w)| *w > 0.0));
        assert!(q.volume_bounds().lo() <= 0.37 && q.volume_bounds().hi() >= 0.37);
        n[axis] = -1.0;
        let other = rule(&Plane { normal: n, offset: -0.37 }, 1);
        assert!((q.volume() + other.volume() - 1.0).abs() < 1e-9);
    }
}
#[test]
fn roots_on_subdivision_planes_do_not_become_ambiguous_or_empty() {
    let q = rule(&Plane { normal: [1.0,0.0,0.0], offset: 0.5 }, 2);
    assert!((q.volume() - 0.5).abs() < 1e-9);
    let full = rule(&Plane { normal: [1.0,0.0,0.0], offset: 2.0 }, 0);
    assert!((full.volume() - 1.0).abs() < 1e-14);
    let empty = rule(&Plane { normal: [1.0,0.0,0.0], offset: -1.0 }, 0);
    assert!(empty.bulk().is_empty());
    assert_eq!(empty.volume().to_bits(), 0.0_f64.to_bits());
}
#[test]
fn curved_domain_volume_is_numerical_and_has_separate_conservative_bounds() {
    let q = rule(&Sphere, 3);
    let exact = (4.0/3.0) * std::f64::consts::PI * 0.35*0.35*0.35;
    assert!((q.volume() - exact).abs()/exact < 0.015, "{} vs {exact}", q.volume());
    assert!(q.volume_bounds().lo() < exact && q.volume_bounds().hi() > exact);
    assert!(q.volume_bounds().hi() - q.volume_bounds().lo() > 1e-6);
    assert!(q.cut_boxes() > 0);
}
#[test]
fn hidden_curved_domain_refuses_when_no_monotone_direction_is_proven() {
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { depth: 0, ..Default::default() }, &mut callback).unwrap();
    assert!(matches!(cut_cell_rules3(&Sphere, unit(), &mut control), Err(QuadratureError3::UnresolvedCell(_))));
}
#[test]
fn budgets_are_cumulative_and_cancellation_does_not_publish_partial_rules() {
    let p = Plane { normal: [1.0,0.0,0.0], offset: 2.0 };
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { max_points: 27, ..Default::default() }, &mut callback).unwrap();
    cut_cell_rules3(&p, unit(), &mut control).unwrap();
    assert!(matches!(cut_cell_rules3(&p, unit(), &mut control), Err(QuadratureError3::PointBudget)));
    assert_eq!(control.work().points, 27);
    let mut callback = |w: fs_cutfem::quad3::QuadratureWork3| {
        if w.field_evaluations > 5 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    };
    let mut control = QuadratureControl3::new(QuadratureOptions3::default(), &mut callback).unwrap();
    assert!(matches!(cut_cell_rules3(&Sphere, unit(), &mut control), Err(QuadratureError3::Cancelled)));
    assert_eq!(control.work().field_evaluations, 6);
}
#[test]
fn deterministic_replay_and_anisotropic_box_measure() {
    let p = Plane { normal: [1.0,0.0,0.0], offset: 0.37 };
    assert_eq!(rule(&p, 2).bulk(), rule(&p, 2).bulk());
    let cell = HexCell::try_new([-1.0,2.0,3.0], [2.0,4.0,7.0]).unwrap();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = QuadratureControl3::new(QuadratureOptions3 { depth: 0, ..Default::default() }, &mut callback).unwrap();
    let q = cut_cell_rules3(&p, cell, &mut control).unwrap();
    assert!((q.volume() - 1.37*2.0*4.0).abs() < 1e-8);
}
