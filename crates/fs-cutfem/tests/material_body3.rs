//! Native density-dependent body loads and their exact occupancy transpose.
use std::{cell::Cell, ops::ControlFlow};

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityError3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solver::LinearOp;

struct Top;
impl CutSdf3 for Top {
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
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn field(cell: HexCell) -> Result<ElasticAcousticState, ElasticityError3> {
    let id = if cell.lo()[0] < 0.5 { 2 } else { 8 };
    let c = std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    Ok(ElasticAcousticState::try_new_mandel(c, f64::from(id), ContentHash([id; 32]), 0.1).unwrap())
}
fn quadrature() -> QuadratureOptions3 {
    QuadratureOptions3 {
        depth: 0,
        root_relative_tolerance: 1e-12,
        ..Default::default()
    }
}
fn cartesian(native: bool) -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(quadrature(), &mut callback).unwrap();
    if native {
        CutElasticity3::build_with_material_field(
            domain(),
            [2; 3],
            &Top,
            ElasticMaterialField3::ByCell(&field),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut control,
        )
        .unwrap()
    } else {
        CutElasticity3::build(
            domain(),
            [2; 3],
            &Top,
            &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(),
            &|p| p[0] == 0.0,
            Default::default(),
            &mut control,
        )
        .unwrap()
    }
}
fn adaptive(mixed: bool) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = if mixed {
        tree.refined(&[*tree.leaves().iter().next().unwrap()], poll)
            .unwrap()
    } else {
        tree
    };
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(quadrature(), &mut callback).unwrap();
    AdaptiveElasticity3::build_with_material_field(
        domain(),
        &tree,
        &Top,
        ElasticMaterialField3::ByCell(&field),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut control,
    )
    .unwrap()
}
fn acceleration(p: [f64; 3]) -> [f64; 3] {
    [0.25 + 0.5 * p[0], -1.0, 2.0 - p[1]]
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn native_density_occupancy_virtual_work_and_transpose_share_one_prepared_load() {
    let mut op = cartesian(true);
    let calls = Cell::new(0);
    let counted = |p| {
        calls.set(calls.get() + 1);
        acceleration(p)
    };
    let stencil = op
        .prepare_material_body_load(&counted, 10_000, poll)
        .unwrap();
    let prepared_calls = calls.get();
    assert!(prepared_calls > 0);
    assert_eq!((stencil.n(), stencil.cells()), (op.n(), op.cells()));
    op.validate_material_body_load(&stencil).unwrap();
    for cell in 0..op.cells() {
        assert_eq!(
            stencil.cell_material_state(cell),
            op.cell_material_state(cell)
        );
    }
    let fractions: Vec<_> = (0..op.cells()).map(|i| 0.2 + 0.07 * i as f64).collect();
    let rhs = stencil.load(&fractions, poll).unwrap();
    let w: Vec<_> = op.nodes().iter().flat_map(|p| [p[0], 0.0, 0.0]).collect();
    let exact: f64 = op
        .cell_keys()
        .iter()
        .zip(&fractions)
        .map(|(key, fraction)| {
            let lo = key[0] as f64 / 2.0;
            let hi = lo + 0.5;
            let average = 0.25 * (lo + hi) / 2.0 + 0.5 * (lo * lo + lo * hi + hi * hi) / 3.0;
            let volume = 0.25 * if key[2] == 0 { 0.5 } else { 0.23 };
            let density = if key[0] == 0 { 2.0 } else { 8.0 };
            fraction * density * volume * average
        })
        .sum();
    close(dot(&w, &rhs), exact, 2e-11);
    for (node, fixed) in op.fixed().iter().enumerate() {
        if *fixed {
            assert_eq!(&rhs[3 * node..3 * node + 3], &[0.0; 3]);
        }
    }
    let dual: Vec<_> = (0..op.n()).map(|i| 0.1 * ((i % 7) as f64 - 3.0)).collect();
    let pullback = stencil.occupancy_pullback(&dual, poll).unwrap();
    close(dot(&dual, &rhs), dot(&fractions, &pullback), 1e-13);
    let direction: Vec<_> = (0..op.cells())
        .map(|i| 0.02 * ((i % 3) as f64 - 1.0))
        .collect();
    let h = 1e-5;
    let plus: Vec<_> = fractions
        .iter()
        .zip(&direction)
        .map(|(r, d)| r + h * d)
        .collect();
    let minus: Vec<_> = fractions
        .iter()
        .zip(&direction)
        .map(|(r, d)| r - h * d)
        .collect();
    let fd = (dot(&dual, &stencil.load(&plus, poll).unwrap())
        - dot(&dual, &stencil.load(&minus, poll).unwrap()))
        / (2.0 * h);
    close(fd, dot(&direction, &pullback), 2e-10);
    op.set_scales(&vec![0.03; op.cells()]).unwrap();
    assert_eq!(stencil.load(&fractions, poll).unwrap(), rhs);
    assert_eq!(
        stencil.load(&vec![0.0; op.cells()], poll).unwrap(),
        vec![0.0; op.n()]
    );
    assert_eq!(calls.get(), prepared_calls);
}

#[test]
fn adaptive_density_load_is_condensed_once_and_remains_bound_after_wrapping() {
    for mixed in [false, true] {
        let op = adaptive(mixed);
        if mixed {
            assert!(op.physical_nodes().len() > op.nodes().len());
        }
        let stencil = op
            .prepare_material_body_load(&acceleration, 100_000, poll)
            .unwrap();
        let rhs = stencil.load(&vec![1.0; op.cells()], poll).unwrap();
        let reference = op
            .body_load(
                &|p| {
                    let density = if p[0] < 0.5 { 2.0 } else { 8.0 };
                    acceleration(p).map(|a| density * a)
                },
                poll,
            )
            .unwrap();
        for (a, b) in rhs.iter().zip(&reference) {
            close(*a, *b, 2e-12);
        }
        let w: Vec<_> = op.nodes().iter().flat_map(|p| [p[0], 0.0, 0.0]).collect();
        // Integral rho*x*(.25+.5*x) dV over the two material regions.
        close(
            dot(&w, &rhs),
            0.73 * (2.0 * 5.0 / 96.0 + 8.0 * 23.0 / 96.0),
            2e-11,
        );
        let fractions: Vec<_> = (0..op.cells())
            .map(|i| 0.3 + 0.02 * (i % 9) as f64)
            .collect();
        let dual: Vec<_> = (0..op.n())
            .map(|i| ((i % 11) as f64 - 5.0) / 11.0)
            .collect();
        close(
            dot(&dual, &stencil.load(&fractions, poll).unwrap()),
            dot(
                &fractions,
                &stencil.occupancy_pullback(&dual, poll).unwrap(),
            ),
            2e-12,
        );
        let foreign = adaptive(mixed);
        assert!(foreign.validate_material_body_load(&stencil).is_err());
        let wrapper = AdaptiveSolveSpace3::jacobi(op, 100_000_000);
        wrapper
            .elasticity()
            .validate_material_body_load(&stencil)
            .unwrap();
        assert_eq!(
            stencil.load(&vec![1.0; stencil.cells()], poll).unwrap(),
            rhs
        );
    }
}

#[test]
fn native_load_preparation_and_actions_refuse_unavailable_or_partial_results() {
    let legacy = cartesian(false);
    let calls = Cell::new(0);
    let counted = |p| {
        calls.set(calls.get() + 1);
        acceleration(p)
    };
    assert!(
        legacy
            .prepare_material_body_load(&counted, 10_000, poll)
            .is_err()
    );
    assert_eq!(calls.get(), 0);
    let op = cartesian(true);
    let scales = op.scales().to_vec();
    assert!(matches!(
        op.prepare_material_body_load(&counted, 10_000, || {
            if calls.get() > 2 {
                ControlFlow::Break(())
            } else {
                poll()
            }
        }),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(
        op.prepare_material_body_load(&acceleration, 0, poll)
            .is_err()
    );
    assert!(
        op.prepare_material_body_load(&|_| [f64::NAN; 3], 10_000, poll)
            .is_err()
    );
    let zero = op
        .prepare_material_body_load(&|_| [0.0; 3], 0, poll)
        .unwrap();
    assert_eq!(
        zero.load(&vec![1.0; op.cells()], poll).unwrap(),
        vec![0.0; op.n()]
    );
    let stencil = op
        .prepare_material_body_load(&acceleration, 10_000, poll)
        .unwrap();
    let fractions = vec![0.4; op.cells()];
    let baseline = stencil.load(&fractions, poll).unwrap();
    for value in [-0.1, 1.1, f64::NAN] {
        let mut invalid = fractions.clone();
        invalid[0] = value;
        assert!(stencil.load(&invalid, poll).is_err());
    }
    assert!(stencil.load(&[], poll).is_err());
    assert!(stencil.occupancy_pullback(&[], poll).is_err());
    assert!(
        stencil
            .occupancy_pullback(&vec![f64::NAN; op.n()], poll)
            .is_err()
    );
    assert!(matches!(
        stencil.load(&fractions, || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(matches!(
        stencil.occupancy_pullback(&baseline, || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(
        cartesian(true)
            .validate_material_body_load(&stencil)
            .is_err()
    );
    assert!(
        adaptive(false)
            .validate_material_body_load(&stencil)
            .is_err()
    );
    let moved = op;
    moved.validate_material_body_load(&stencil).unwrap();
    assert_eq!(moved.scales(), scales);
    assert_eq!(stencil.load(&fractions, poll).unwrap(), baseline);
}

#[test]
fn body_fingerprint_binds_the_discrete_map_and_sources_across_rebuilds() {
    let first = cartesian(true);
    let rebuilt = cartesian(true);
    let load = first
        .prepare_material_body_load(&acceleration, 10_000, poll)
        .unwrap();
    let same = rebuilt
        .prepare_material_body_load(&acceleration, 10_000, poll)
        .unwrap();
    assert_eq!(load.fingerprint(), same.fingerprint());
    assert!(
        rebuilt.validate_material_body_load(&load).is_err(),
        "content equality must not weaken the native owner binding"
    );
    let changed = rebuilt
        .prepare_material_body_load(&|p| acceleration(p).map(|a| 2.0 * a), 10_000, poll)
        .unwrap();
    assert_ne!(load.fingerprint(), changed.fingerprint());

    // With zero acceleration the complete sparse maps agree, so only source
    // attribution can distinguish these otherwise identical prepared values.
    let relabelled_field = |cell| {
        let source = field(cell)?;
        Ok(ElasticAcousticState::try_new_mandel(
            *source.stiffness_mandel_pa(),
            source.density_kg_m3(),
            ContentHash([99; 32]),
            source.maximum_absolute_mandel_strain(),
        )
        .unwrap())
    };
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(quadrature(), &mut callback).unwrap();
    let relabelled = CutElasticity3::build_with_material_field(
        domain(),
        [2; 3],
        &Top,
        ElasticMaterialField3::ByCell(&relabelled_field),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut control,
    )
    .unwrap();
    let zero = first
        .prepare_material_body_load(&|_| [0.0; 3], 0, poll)
        .unwrap();
    let other_zero = relabelled
        .prepare_material_body_load(&|_| [0.0; 3], 0, poll)
        .unwrap();
    assert_eq!(
        zero.load(&vec![1.0; first.cells()], poll).unwrap(),
        other_zero
            .load(&vec![1.0; relabelled.cells()], poll)
            .unwrap()
    );
    assert_ne!(zero.fingerprint(), other_zero.fingerprint());
}
