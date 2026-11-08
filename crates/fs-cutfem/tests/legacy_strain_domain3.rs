//! Original isotropic strain ceilings must survive CutFEM preparation and solves.
use std::ops::ControlFlow;

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticityError3, ElasticityFormulation3, ElasticityOptions3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;

struct Full;
impl CutSdf3 for Full {
    fn value(&self, _: [f64; 3]) -> f64 {
        -1.0
    }
    fn enclose(&self, _: [f64; 3], _: [f64; 3]) -> Interval {
        Interval::new(-1.0, -1.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], _: HeightAxis) -> Interval {
        Interval::new(0.0, 0.0)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn options(formulation: ElasticityFormulation3) -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation,
        ghost_gamma: 0.0,
        ..Default::default()
    }
}
fn cartesian(limit: f64, formulation: ElasticityFormulation3, fix_face: bool) -> CutElasticity3 {
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(Default::default(), &mut checkpoint).unwrap();
    CutElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [1; 3],
        &Full,
        &IsotropicElastic::new(10.0, 0.0, limit).unwrap(),
        &|p| if fix_face { p[0] == 0.0 } else { p == [0.0; 3] },
        options(formulation),
        &mut control,
    )
    .unwrap()
}
fn adaptive(limit: f64, formulation: ElasticityFormulation3) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 3, 128).unwrap();
    let mark = *tree
        .leaves()
        .iter()
        .find(|c| c.index() == [1, 0, 0])
        .unwrap();
    let tree = tree.refined(&[mark], poll).unwrap();
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Full,
        &IsotropicElastic::new(10.0, 0.0, limit).unwrap(),
        &|p| p[0] == 0.0,
        options(formulation),
        &mut control,
    )
    .unwrap()
}
fn field(nodes: &[[f64; 3]], f: impl Fn([f64; 3]) -> [f64; 3]) -> Vec<f64> {
    nodes.iter().flat_map(|p| f(*p)).collect()
}

#[test]
fn isotropic_limit_is_the_full_tensor_norm_with_both_shear_entries() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let mut op = cartesian(0.1, formulation, false);
        let accepted = field(op.nodes(), |p| [0.06 * p[0], 0.06 * p[1], 0.0]);
        let excessive = field(op.nodes(), |p| [0.08 * p[0], 0.08 * p[1], 0.0]);
        let shear = field(op.nodes(), |p| [0.08 * p[1], 0.08 * p[0], 0.0]);
        op.admit_material_domain(&accepted, poll).unwrap();
        let points = op.bulk_stress(&accepted, 10_000, poll).unwrap();
        let observations = op.bulk_strain_stress(&accepted, 10_000, poll).unwrap();
        for (observation, point) in observations.iter().zip(&points) {
            assert_eq!(&observation.response, point);
        }
        let cotangents = vec![[1.0; 6]; points.len()];
        for u in [&excessive, &shear] {
            assert!(op.admit_material_domain(u, poll).is_err());
            assert!(op.bulk_stress(u, 10_000, poll).is_err());
            assert!(op.bulk_strain_stress(u, 10_000, poll).is_err());
            assert!(op.bulk_stress_pullback(u, &cotangents, poll).is_err());
        }
        // A reference transpose is an algebraic action without a physical primal.
        op.reference_bulk_stress_pullback(&cotangents, poll)
            .unwrap();
        // A stiffness interpolation does not change the physical strain domain.
        op.set_scales(&vec![0.001; op.cells()]).unwrap();
        assert!(op.admit_material_domain(&excessive, poll).is_err());
    }
}

#[test]
fn isotropic_domain_covers_corner_peaks_and_correlated_load_combinations() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = cartesian(0.1, formulation, false);
        let corner_peak = field(op.nodes(), |p| [0.085 * p[0] * p[1], 0.0, 0.0]);
        // Interior samples are admissible, but epsilon_xx and epsilon_xy
        // together exceed the Frobenius ceiling near the corner.
        op.bulk_stress(&corner_peak, 10_000, poll).unwrap();
        assert!(op.admit_material_domain(&corner_peak, poll).is_err());

        let x = field(op.nodes(), |p| [0.09 * p[0], 0.0, 0.0]);
        let y = field(op.nodes(), |p| [0.0, 0.09 * p[1], 0.0]);
        op.admit_load_space_material_domain(&[&x, &y], poll)
            .unwrap();
        assert!(
            op.admit_load_space_material_domain(&[&x, &x], poll)
                .is_err()
        );
        let zero = vec![0.0; op.n()];
        assert!(
            op.affine_bulk_stress(&zero, &[&x, &x], 20_000, poll)
                .is_err()
        );
        let offset = field(op.nodes(), |p| [0.06 * p[0], 0.0, 0.0]);
        let direction = field(op.nodes(), |p| [0.0, 0.07 * p[1], 0.0]);
        op.admit_affine_material_domain(&offset, &[&direction], poll)
            .unwrap();
        op.affine_bulk_stress(&offset, &[&direction], 20_000, poll)
            .unwrap();
    }
}

#[test]
fn adaptive_material_domain_uses_the_retained_limit_after_hanging_reconstruction() {
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = adaptive(0.1, formulation);
        assert!(op.physical_nodes().len() > op.nodes().len());
        let accepted = field(op.nodes(), |p| [0.09 * p[0], 0.0, 0.0]);
        let excessive = field(op.nodes(), |p| [0.11 * p[0], 0.0, 0.0]);
        op.admit_material_domain(&accepted, poll).unwrap();
        op.bulk_stress(&accepted, 100_000, poll).unwrap();
        op.bulk_strain_stress(&accepted, 100_000, poll).unwrap();
        assert!(op.admit_material_domain(&excessive, poll).is_err());
        assert!(op.bulk_stress(&excessive, 100_000, poll).is_err());
        assert!(op.bulk_strain_stress(&excessive, 100_000, poll).is_err());
        assert!(
            op.admit_load_space_material_domain(&[&accepted, &accepted], poll)
                .is_err()
        );
    }
}

macro_rules! physical_solve_checks {
    ($op:expr) => {{
        let op = $op;
        for (strain, admissible) in [(0.005, true), (0.02, false)] {
            let expected = field(op.nodes(), |p| [strain * p[0], 0.0, 0.0]);
            let mut force = vec![0.0; op.n()];
            op.apply(&expected, &mut force);
            let algebraic = op
                .solve_controlled(&force, 1e-11, 20_000, 16, |_| poll())
                .unwrap();
            let physical = op.solve_physical_controlled(&force, 1e-11, 20_000, 16, |_| poll());
            if admissible {
                assert_eq!(physical.unwrap().coefficients(), algebraic.coefficients());
            } else {
                assert!(matches!(physical, Err(ElasticityError3::Invalid(_))));
            }
            assert!(matches!(
                op.solve_physical_controlled(&force, 1e-11, 20_000, 16, |_| ControlFlow::Break(())),
                Err(ElasticityError3::Cancelled)
            ));
        }
    }};
}

#[test]
fn physical_solves_publish_only_admitted_fields_and_keep_adjoint_solves_available() {
    physical_solve_checks!(cartesian(0.01, ElasticityFormulation3::Displacement, true));
    physical_solve_checks!(adaptive(0.01, ElasticityFormulation3::Displacement));
}

#[test]
fn isotropic_admission_removes_rigid_translation_before_differentiating() {
    struct Half;
    impl CutSdf3 for Half {
        fn value(&self, p: [f64; 3]) -> f64 {
            p[2] - 0.73
        }
        fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
            Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
        }
        fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
            let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
            Interval::new(d, d)
        }
    }
    let mut checkpoint = |_| poll();
    let mut control = QuadratureControl3::new(Default::default(), &mut checkpoint).unwrap();
    let op = CutElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [0.3, 0.7, 0.9]).unwrap(),
        [2; 3],
        &Half,
        &IsotropicElastic::new(10.0, 0.0, 1e-10).unwrap(),
        &|_| false,
        &|_, _| true,
        Default::default(),
        Default::default(),
        Default::default(),
        &mut control,
    )
    .unwrap();
    let translated = field(op.nodes(), |_| [1e12, -2e12, 3e12]);
    op.admit_material_domain(&translated, poll).unwrap();
    let observations = op.bulk_strain_stress(&translated, 100_000, poll).unwrap();
    assert!(observations.iter().all(|point| point.strain == [0.0; 6]
        && point.constitutive_strain == [0.0; 6]
        && point.response.stress == [0.0; 6]
        && point.energy_density == 0.0));
}
